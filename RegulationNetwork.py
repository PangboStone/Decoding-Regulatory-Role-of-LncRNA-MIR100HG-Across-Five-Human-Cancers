import torch
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from torch_geometric.data import Data, DataLoader
from torch_geometric.nn import GATConv, global_mean_pool
from node2vec import Node2Vec
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


def load_data():

    gene_expr = pd.read_csv("Gene_expression.csv", index_col=0)  # Genes x Samples
    corr_matrix = pd.read_csv("correlations.csv", index_col=0)  # TFs x MIR100HG
    

    high_corr_tfs = corr_matrix[abs(corr_matrix['MIR100HG']) > 0.5].index.tolist()
    
    # Building global network
    G = nx.Graph()
    G.add_node("MIR100HG", features=gene_expr.loc["MIR100HG"].values)
    
    for tf in high_corr_tfs:
        G.add_node(tf, features=gene_expr.loc[tf].values)
        G.add_edge("MIR100HG", tf, weight=corr_matrix.loc[tf, 'MIR100HG'])
    
    return G

# Generating Subgraphs
def generate_subgraphs(G, num_samples=1000, noise_ratio=0.2):
    subgraphs = []
    nodes = list(G.nodes())
    
    for _ in range(num_samples):
        
        tf = np.random.choice([n for n in nodes if n != "MIR100HG"])
        
        # Constructing initial subgraph
        base_subgraph = nx.ego_graph(G, tf, radius=2)
        
        # Add Noise Edge
        noise_edges = []
        for _ in range(int(noise_ratio*base_subgraph.number_of_edges())):
            n1, n2 = np.random.choice(nodes, 2, replace=False)
            noise_edges.append((n1, n2))
            
        noisy_subgraph = nx.Graph(base_subgraph)
        noisy_subgraph.add_edges_from(noise_edges)
        
        # Based on initial correlation generate Label
        label = 1 if G.edges[("MIR100HG", tf)]['weight'] > 0 else 0
        subgraphs.append((noisy_subgraph, label))
    
    return subgraphs


class FeatureEngineer:
    def __init__(self, G):
        self.G = G
        self.node2vec = Node2Vec(G, dimensions=64, walk_length=30, num_walks=200)
        
    def get_features(self, subgraph):
        # 结构嵌入
        model = self.node2vec.fit()
        struct_emb = {node: model.wv[node] for node in subgraph.nodes()}
        
        # 节点特征矩阵
        features = []
        for node in subgraph.nodes():
            expr_feat = self.G.nodes[node]['features']
            emb_feat = struct_emb[node]
            combined = np.concatenate([expr_feat, emb_feat])
            features.append(combined)
        
        return np.array(features)


class RegulatoryGNN(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim=128):
        super().__init__()
        self.conv1 = GATConv(input_dim, hidden_dim)
        self.conv2 = GATConv(hidden_dim, hidden_dim)
        self.classifier = torch.nn.Linear(hidden_dim, 1)
        
    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index).relu()
        x = global_mean_pool(x, data.batch)
        return torch.sigmoid(self.classifier(x))


class EnsembleRegulatoryNetwork:
    def __init__(self, num_models=5):
        self.models = [RegulatoryGNN(input_dim=100) for _ in range(num_models)]  # 假设特征维度100
        self.ensemble_weights = None
        
    def train_single(self, model, loader, epochs=100):
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        criterion = torch.nn.BCELoss()
        
        for epoch in range(epochs):
            for data in loader:
                optimizer.zero_grad()
                out = model(data)
                loss = criterion(out, data.y)
                loss.backward()
                optimizer.step()
                
    def train_ensemble(self, subgraphs):
        pyg_data = []
        fe = FeatureEngineer(G)
        
        for sg, label in subgraphs:
            features = fe.get_features(sg)
            edge_index = torch.tensor(list(sg.edges())).t().contiguous()
            data = Data(x=torch.FloatTensor(features),
                        edge_index=edge_index,
                        y=torch.FloatTensor([label]))
            pyg_data.append(data)
            
        train_data, val_data = train_test_split(pyg_data, test_size=0.2)
        loader = DataLoader(train_data, batch_size=32, shuffle=True)
        
        # train
        for i, model in enumerate(self.models):
            print(f"Training model {i+1}")
            self.train_single(model, loader)
            
        # Learning weights integration
        X_meta, y_meta = [], []
        for data in val_data:
            preds = [m(data).item() for m in self.models]
            X_meta.append(preds)
            y_meta.append(data.y.item())
            
        self.ensemble_weights = LogisticRegression().fit(X_meta, y_meta)
        
    def predict(self, subgraph):
        with torch.no_grad():
            preds = [m(subgraph).item() for m in self.models]
            return self.ensemble_weights.predict_proba([preds])[0][1]

#  Network construction and visualisation
def build_final_network(G, ensemble_model, threshold=0.7):
    final_net = nx.Graph()
    fe = FeatureEngineer(G)
    
    for node in G.nodes():
        final_net.add_node(node)
        
    for u in G.nodes():
        if u == "MIR100HG":
            continue
        subgraph = nx.ego_graph(G, u, radius=2)
        prob = ensemble_model.predict(subgraph)
        if prob > threshold:
            final_net.add_edge("MIR100HG", u, weight=prob)
    
    # visualization
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(final_net)
    nx.draw(final_net, pos, with_labels=True, 
            node_color=['red' if n=='MIR100HG' else 'skyblue' for n in final_net.nodes()],
            edge_color=[final_net[u][v]['weight'] for u,v in final_net.edges()],
            width=2, edge_cmap=plt.cm.Blues)
    plt.title("MIR100HG Regulatory Network")
    plt.show()
    
    return final_net

if __name__ == "__main__":

    G = load_data()
    

    subgraphs = generate_subgraphs(G)
    

    ensemble = EnsembleRegulatoryNetwork(num_models=5)
    ensemble.train_ensemble(subgraphs)
    

    final_net = build_final_network(G, ensemble)
    

    nx.write_adjlist(final_net, "MIR100HG_regulatory_network.adjlist")