#include <iostream>
#include <queue>
#include <vector>
using namespace std;
typedef vector<int> valNode;
typedef vector<valNode> adyacencias;

int PrimsMST(int sourceNode, vector<adyacencias>& graph, int K)
{
    //Guardar detalles del nodo.
    priority_queue<valNode, vector<valNode>, greater<valNode>> k;
    int count = 0;
    vector<int> aux = { 0,sourceNode };
    k.push(aux);
    bool* nodesAdded = new bool[graph.size()];
    memset(nodesAdded, false, sizeof(bool) * graph.size());
    int mst_tree_cost = 0;

    while (count != K)
    {
        // Nodo con mínimo costo
        valNode itemNode;
        itemNode = k.top();
        k.pop();
        int Node = itemNode[1];
        int Cost = itemNode[0];

        //Checar si el nodo ya se añadió
        if (!nodesAdded[Node]) 
        {
            mst_tree_cost += Cost;
            count++;
            if (count == K)
                break;
            nodesAdded[Node] = true;

            //Nodos vecinos quitados de priority queque
            for (auto& node_cost : graph[Node]) 
            {
                int adjacency_node = node_cost[1];
                if (nodesAdded[adjacency_node] == false) 
                {
                    k.push(node_cost);
                }
            }
        }
    }
    delete[] nodesAdded;
    return mst_tree_cost;
}


int main() 
{
    adyacencias fromNode_0_in_graph_1 = { {1,1}, {2,2}, {1,3}, {1,4}, {2,5}, {1,6} };
    adyacencias fromNode_1_in_graph_1 = { {1,0}, {2,2}, {2,6} };
    adyacencias fromNode_2_in_graph_1 = { {2,0}, {2,1}, {1,3} };
    adyacencias fromNode_3_in_graph_1 = { {1,0}, {1,2}, {2,4} };
    adyacencias fromNode_4_in_graph_1 = { {1,0}, {2,3}, {2,5} };
    adyacencias fromNode_5_in_graph_1 = { {2,0}, {2,4}, {1,6} };
    adyacencias fromNode_6_in_graph_1 = { {1,0}, {2,2}, {1,5} };

    int num_of_nodes = 7; // Total Nodes (0 to 6)
    vector<adyacencias> primsgraph;
    primsgraph.resize(num_of_nodes);
    primsgraph[0] = fromNode_0_in_graph_1;
    primsgraph[1] = fromNode_1_in_graph_1;
    primsgraph[2] = fromNode_2_in_graph_1;
    primsgraph[3] = fromNode_3_in_graph_1;
    primsgraph[4] = fromNode_4_in_graph_1;
    primsgraph[5] = fromNode_5_in_graph_1;
    primsgraph[6] = fromNode_6_in_graph_1;

    // As we already know, we have to choose the source vertex,
    // so we start from the vertex 0 node.
    cout << "k-mst : " "" << PrimsMST(3, primsgraph, 3) << std::endl;
    return 0;
}
