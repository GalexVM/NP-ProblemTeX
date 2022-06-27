#include <iostream>
#include <unordered_set>
#include <set>
#include <list>
#include <queue>
#include <bitset> 
#include "grafos.h"
#include "abstract.h"
using namespace std;

typedef Graf<int, int> G;
//typedef unordered_set<Arista<G>> HashSet;

typedef bitset<4> BitSet;

void PrintMatrix(vector<vector<int>> mat) {
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			cout << mat[i][j] << "\t";
		}
		cout << endl;
	}
}

//////////////FUERZA BRUTA: PRIM//////////////////
int minKey(vector<int> key, vector<bool>mstSet)
{
	// Initialize min value
	int min = INT_MAX, min_index;

	for (int v = 0; v < key.size(); v++)
		if (mstSet[v] == false && key[v] < min)
			min = key[v], min_index = v;

	return min_index;
}

// A utility function to print the
// constructed MST stored in parent[]
void printMST(vector<int> parent, vector<vector<int>> graph)
{
	cout << "Arista \tPeso\n";
	for (int i = 1; i < graph.size(); i++)
		cout << parent[i] << " - " << i << " \t" << graph[i][parent[i]] << " \n";
}

// Function to construct and print MST for
// a graph represented using adjacency
// matrix representation
void primMST(vector<vector<int>> graph)
{
	int V = graph.size();
	// Array to store constructed MST
	vector<int> parent(V);

	// Key values used to pick minimum weight edge in cut
	vector<int> key(V);

	// To represent set of vertices included in MST
	vector<bool> mstSet(V);

	// Initialize all keys as INFINITE
	for (int i = 0; i < V; i++)
		key[i] = INT_MAX, mstSet[i] = false;

	// Always include first 1st vertex in MST.
	// Make key 0 so that this vertex is picked as first vertex.
	key[0] = 0;
	parent[0] = -1; // First node is always root of MST

	// The MST will have V vertices
	for (int count = 0; count < V - 1; count++)
	{
		// Pick the minimum key vertex from the
		// set of vertices not yet included in MST
		int u = minKey(key, mstSet);

		// Add the picked vertex to the MST Set
		mstSet[u] = true;

		// Update key value and parent index of
		// the adjacent vertices of the picked vertex.
		// Consider only those vertices which are not
		// yet included in MST
		for (int v = 0; v < V; v++)

			// graph[u][v] is non zero only for adjacent vertices of m
			// mstSet[v] is false for vertices not yet included in MST
			// Update the key only if graph[u][v] is smaller than key[v]
			if (graph[u][v] && mstSet[v] == false && graph[u][v] < key[v])
				parent[v] = u, key[v] = graph[u][v];
	}

	printMST(parent, graph);
}





//////////////APROXIMADO: PRIM OPTIMIZADO + BRANCH & BOUND//////////////////



//KMST
class KMST {
private:
	vector<vector<Arista<G>>*> edgesFromNode;

	//vector<vector<Arista<G>>> edgesFromNode;
  //  Arista<G>* edgesFromNode;
	  //unordered_set<HashSet> visited;
	set<HashSet> visited;
	// HashSet visited;
	vector<int> minSum;
	int numNodes;
	int numEdges;
	int k = 0;
	int minWeight = INT_MAX;
	int kEdges;
	int limit;
	int abort;
	HashSet edges;
	bool limited;

public:
	//------------------------------------------------------------------
	bool hasNoCircle(BitSet used, int node1, int node2) {
		if (used[node1] && used[node2]) {
			return false;
		}
		return true;
	}
	//-----------------------------------------
	KMST(int numNodess, int numEdgess, HashSet edges, int k) {
		this->numNodes = numNodess;
		this->numEdges = numEdgess;
		this->edges = edges;
		this->k = k;
		this->kEdges = k - 1;
		// visited = new set<HashSet>;
		edgesFromNode.resize(numNodes);
		minSum.resize(k);
		this->abort = 0;
		this->limited = true;

		// sane limit for most graphs
		this->limit = 2 * numNodes * numNodes;


		// PriorityQueue for the k cheapest edges	
		priority_queue<Arista<G>>* min = new priority_queue<Arista<G>>;

		// Create data structure (adjacency list)
		for (Arista<G> t : edges) {

			if (edgesFromNode[t.node1] == 0) {
				edgesFromNode[t.node1] = new vector<Arista<G>>(numNodes);
			}
			if (edgesFromNode[t.node2] == 0) {
				edgesFromNode[t.node2] = new vector<Arista<G>>(numNodes);
			}
			edgesFromNode[t.node1]->push_back(t);
			edgesFromNode[t.node2]->push_back(t);
			min->push(t);
		}

		// i use the sum of the k - |V| cheapest edges to determine if a given
		// graph could ever be better than the minWeight
		minSum[0] = 0;
		for (int i = 1; i < k; i++) {
			minSum[i] = minSum[i - 1] + min->top().m_v;
			min->pop();
		}
	}

	//---------------------------------------------------------------
	void run() {
		constructMST();

	}

	//------------------------------------------------------------------------------
	void constructMST() {
    vector<Arista<G>> aux(numNodes);
		priority_queue<Arista<G>> q(aux.begin(),aux.end());
		int t;

		// adds the sums of all nodes to the PriorityQueue in reverse order
		for (int i = 0; i < numNodes; i++) {
			t = getBestEdge(i);
			if (t != INT_MAX) {
				Arista<G> a(t);
				a.node1 = i;
				a.node2 = -1;
        Arista<G> ari(i, -1, t);
				q.push(ari);
			}
		}

		// clone node queues
		priority_queue<Arista<G>> q_prim = q;
		priority_queue<Arista<G>> q_limited = q;

		// fast upper bound
		// estimate to determine a good upper bound - modified prims algorithm
		// no backtracking
		while (!q_prim.empty()) {
      vector<Arista<G>> aux(k);
      vector<Arista<G>> aux2(numEdges);
			firstEstimate(*new HashSet(aux.begin(), aux.end()), q_prim.top().node2, 0,
					*new priority_queue<Arista<G>>(aux2.begin(), aux2.end()), *new BitSet(numNodes), 0);
			q_prim.pop();
		}

		// limited enum
		// beginning with the best node it enumerates all possible solutions
		// (branch) until the recursion limit is reached and cuts if the graph
		// is useless
		while (!q_limited.empty()) {
      HashSet n;
      //comente esto  
     // priority_queue<Arista<G>> q;
			addNodes(n, q_limited.top().node2, 0, q,
					*new BitSet(numNodes), 0);
			q_limited.pop();
			abort = 0;
		}

		visited.clear();
		limited = false;
		limit = INT_MAX;

		// full enum
		// beginning with the best node it enumerates all possible solutions
		// (branch) and cuts if the graph is useless
		while (!q.empty()) {
      HashSet n;
      priority_queue<Arista<G>> p;
			addNodes(n, q.top().node2, 0, p, *new BitSet(numNodes), 0);
			q.pop();
		}
    cout<<"finish"<<endl;
    
	}




	//---------------------------------------------------
	int getBestArista(int node, vector<vector<int>> mat) {
		int ret = 0;
		for (int e : mat[node]) {
			ret += e;
		}
		return ret * -1;

	}

	//---------------------------------------------------------------
	bool find(priority_queue<Arista<G>> e, const int& val) const
	{
		while (!e.empty()) {
			if (e.top().m_v == val) return true;
			e.pop();
		}
		return false;
	}

	//------------------------------------------------------------------------
	bool find(set<HashSet> e, HashSet adj) const
	{
		for (auto& a : e) {
			if (a == adj)
				return true;
		}
		return false;
	}

	//------------------------------------------------------------------------ falta acomodar*
	void addToQueue(priority_queue<Arista<G>> e, int node, BitSet used, int w,
		int numEdges) {
		// we iterate through all adjacent nodes using the adjacency list
		for (Arista<G> ite : (*edgesFromNode[node])) {
			// if the expaning node == ite.m_nodes[0] we check if nodes[1] is used to
			// prevent a circle and vice verca; the weight of the current graph
			// + the weight of a new edge + the weight of the kEdges - |E|
			// cheapest edges must be < minWeight; we do not allow duplicates in
			// the todo-edge-list
			if (!used[node == ite.node1 ? ite.node2 : ite.node2]
				&& w + ite.m_v + minSum[kEdges - numEdges - 1] < minWeight
				&& !find(e, ite.m_v)) {
				e.push(ite);
			}
		}
	}

	//------------------------------------------------------------------------
	void firstEstimate(HashSet e, int node, int cweight,
		priority_queue<Arista<G>> p, BitSet used, int numAristas) {
		Arista<G> t;
		int w, newNode;
		bool abort = false, wasEmpty, solutionFound;

		// adds elements adj. to node to the edge-queue 
		addToQueue(p, node, used, cweight, numAristas);

		while (!p.empty() && !abort) {
			t = p.top();
			p.pop();

			// if a given node has a higher weight than minWeight we can ignore
			// it entirely - very unlikely
			if (t.m_v >= minWeight) {
				edgesFromNode[t.node1]->erase((edgesFromNode[t.node1]->begin() + t.node2));
				edgesFromNode[t.node2]->erase((edgesFromNode[t.node2]->begin() + t.node1));

			}
			else {
				w = cweight + t.m_v;

				// circle check
				if (hasNoCircle(used, t.node1, t.node2)) {
					// make sure to quit the loop
					abort = true;

					if (used[t.node1]) {
						// m_nodes[0] is already in use => nodes[1] is new
						newNode = t.node2;
						node = t.node1;
					}
					else {
						newNode = t.node1;
						node = t.node2;
					}

					// add edge to solution
					e.insert(t);

					wasEmpty = false;
					solutionFound = false;

					if (used.none()) {
						// first edge
						used.set(newNode);
						used.set(node);
						wasEmpty = true;
					}
					else {
						used.set(newNode);
					}

					int size = used.count();

					// if |V| = k and the solution is better than minWeight, we
					// update our best solution
					if (size == k && w < minWeight) {
						minWeight = w;
						AbstractKMST o; //implementR UPDATEsolution
						o.setSolution(w, e);
					}
					else if (size < k) {
						// we need to add more edges
						firstEstimate(e, newNode, w, p, used, numAristas + 1);
					}
					// removes the used nodes
					if (!solutionFound) {
						used.reset(newNode);
						if (wasEmpty) {
							used.reset(node);
						}
					}
				}
			}
		}
	}

	//------------------------------------------------------------------
	int getBestEdge(int node) {
		int ret = 0;
		for (Arista<G> e : *edgesFromNode[node]) {
			ret += e.m_v;
		}
		return ret * -1;

	}


	//------------------------------------------------------------------
	void addNodes(HashSet e, int node, int cweight, priority_queue<Arista<G>> p, BitSet used, int numAristas) {


		if (limited)
			abort++;

		Arista<G> t;
		vector<Arista<G>> aux(2 * k);
		HashSet* temp = new HashSet(aux.begin(), aux.end());
		int w, newNode, size;
		bool wasEmpty, solutionFound;

		// clone
		if (!p.empty()) {
			p = *new priority_queue<Arista<G>>(p);
		}
		else {
			p = *new priority_queue<Arista<G>>();
		}

		 if (used != NULL) {
		 	used = used;
		 }

		if (!e.empty()) {
			temp->insert(e.begin(), e.end());
		}

		// expand node
		addToQueue(p, node, used, cweight, numAristas);

		while (!p.empty() && abort < limit) {
			t = p.top();
			p.pop();
			w = cweight + t.m_v;

			// if the weight of the current graph plus the weight of the (k -
			// |V|) cheapest edges is greater than minWeight, we can abort. we
			// also stop enumerating if the current graph has been expanded
			// before

			if (w + minSum[kEdges - numAristas - 1] < minWeight
				&& !find(visited, *temp)) {

				// circle check
				if (hasNoCircle(used, t.node1, t.node2)) {
					if (used[t.node1]) {
						// m_nodes[0] is part of the graph => nodes[1] is new
						newNode = t.node2;
						node = t.node1;
					}
					else {
						newNode = t.node1;
						node = t.node2;
					}

					temp->insert(t);

					wasEmpty = false;
					solutionFound = false;
					if (used.none()) {
						// first edge
						used.set(newNode);
						used.set(node);
						wasEmpty = true;
					}
					else {
						used.set(newNode);
					}

					// number of used nodes
					size = used.count();

					if (size == k) {
						// new best solution found
						updateSolution(*temp, w);
						solutionFound = true;
						abort = 0;
					}
					else {
						// we need to expand more
						addNodes(*temp, newNode, w, p, used, numAristas + 1);
						// if the graph contains 2 nodes we save it to prevent
						// the repeated enumeration of the same solutions
						if (size == 2) {
							visited.insert(*temp);
						}

						// revert to starting solution
						vector<Arista<G>> helper(k); //esto solo es para poder inicializar el set con un tamaño especifico
						temp = new HashSet(helper.begin(), helper.end());
						if (!e.empty()) {
							temp->insert(e.begin(), e.end());
						}
					}
					// clear nodes
					if (!solutionFound) {
						used.reset(newNode);
						if (wasEmpty) {
							used.reset(node);
						}
					}
				}
			}
        //si encuentra la solucion
			else {
        cout<<"solucion encontrada"<<endl;
				break;
			}
		}
	}

	//------------------------------------------------------------------
	void updateSolution(HashSet minSet, int min) {
		minWeight = min;
		AbstractKMST a;
		a.setSolution(min, minSet);
		cout<<min<<endl;

	}
	//-----------------------------------------
};


int main() 
{
	G LOL;
	LOL.insertNode(1);
	LOL.insertNode(2);
	LOL.insertNode(3);
	LOL.insertNode(4);
	LOL.insertNode(5);
	LOL.insertNode(6);
	LOL.insertArista(1, 2, 3, false);
	LOL.insertArista(5, 4, 4, false);
	LOL.insertArista(5, 3, 3, false);
	LOL.insertArista(4, 2, 1, false);
  LOL.insertArista(4, 6, 6, false);
	vector<vector<int>> a = AdjacencyFromGraph(LOL);
	PrintMatrix(a);

	HashSet ar = HashFromGraph(LOL);  
  for (auto& ari : ar){
    cout<<ari.m_nodos[0]->id<<" "<<ari.m_v<<endl;
  }
	//KMST(int numNodess, int numEdgess, HashSet edges, int k) 
	// KMST b(6, 4, ar, 2);
	// b.run();

}