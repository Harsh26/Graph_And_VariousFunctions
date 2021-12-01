// Graphs.cpp : Defines the entry point for the console application.


#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <limits.h>

using namespace std;

/* ***************All Structures Start ********************* */

// Adjacency Node Structure
class adjNode 
{
	public:
			int val, cost;
		    adjNode* next;

			//DFS
			char type; // t->tree f->forward b->back c->cross 
};



// Graph edge structure
class graphEdge 
{
	public:
			int start_ver, end_ver, weight;
};

// Directed Graph Structure
class DiaGraph
{
	 int N;  // number of nodes in the graph

    // insert new nodes into adjacency list from given graph
    adjNode* getAdjListNode(int value, int weight, adjNode* head)   
	{
        adjNode* newNode = new adjNode;
        newNode->val = value;
        newNode->cost = weight;
         
        newNode->next = head;   // point new node to current head
        return newNode;
    }

public:

    adjNode **head; //adjacency list as an array of pointers

    // Constructor for initilisation
    DiaGraph(vector<graphEdge> edges, int n, int N)  
	{
        // allocate new node
        head = new adjNode*[N]();

        this->N = N;

        // initialize head pointer for all vertices
        for (int i = 0; i < N; ++i)
            head[i] = NULL;

        // construct directed graph by adding edges
        for (int i = 0; i < n; i++)  
		{
            int start_ver = edges[i].start_ver;
            int end_ver = edges[i].end_ver;
            int weight = edges[i].weight;
            // insert in the beginning
            adjNode* newNode = getAdjListNode(end_ver, weight, head[start_ver]);
             
            // point head pointer to new node
            head[start_ver] = newNode;
          }

    }

      // Destructor
     ~DiaGraph() 
	 {
		for (int i = 0; i < N; i++)
			delete[] head[i];

        delete[] head;
     }
};

// Linked list structure to give Topological Sequence
class TopologicalOrder
{
	public:
		int vertex;
		TopologicalOrder *next;
};



// Min heap node structure
class MinHeapNode
{
	public:
			int  v;
			int dist; // Dist from source act as priority. The lesser the dist, more will be priority
};

// Structure to represent a min heap
class MinHeap
{

	public: 
			// Number of heap nodes present currently
			int size;    
   
			// Capacity of min heap
			int capacity; 
   
			// For decreaseKey()
			int *pos;   


			// Heap array
            MinHeapNode **arr;
};


/* ***************All Structures End ********************* */



/* ***************All Utility Functions Start ********************* */

// Add Edge Utility function
void addEdge(vector<graphEdge> &ge,int u,int v,int wt)
{
	graphEdge e;

	e.start_ver=u;
	e.end_ver=v;
	e.weight=wt;

	ge.push_back(e);
}

void createOwnFile(string filename)
{
	filename.insert(0,"./input/");
	ofstream MyFile(filename);

	int V,E;

	int u,v,w;

	cout<<"\n\nEnter number of vertices: ";
	cin>>V;

	cout<<"\n\nEnter number of edges: ";
	cin>>E;

	MyFile << V << " " << E << "\n";

	if(V==0 || V==1)
	{
		cout<<"\nLow number of vertices. No edge possible\n";
		return;
	}

	cout<<"\n\nPlease enter data with utmost care otherwise program may not run.\n\n";
	cout<<"\n\nVertices are numbered from "<<0<<" - "<<V-1<<" both inclusive\n\n";
	cout<<"\n\nWhen entering vertex number it should be between above mentioned limits only\n\n";

	for(int i=0;i<E;i++)
	{
		cout<<"\n\nEdge "<<i+1<<" details:\n\n";

		cout<<"Enter vertex from which edge would start [Vertex is an integer from "<< 0 <<" - "<<V-1<<"]: ";
		cin>>u;

		cout<<"Enter vertex to which edge would end [Vertex is an integer from "<< 0 <<" - "<<V-1<<"]: ";
		cin>>v;

		if(u<=-1  || u>V-1 || v<=-1 || v>V-1)
		{
			cout<<"\n\nYou entered something wrong. Enter again\n\n";
			i--;
			continue;
		}

		cout<<"Enter Weight of edge: ";
		cin>>w;

		MyFile<< u <<" "<<v<< " "<<w<<"\n";

		cout<<endl;
	}

	cout<<"\n\nYour file with name & location "<< filename<<" created Successfully!!!\n\n";

	// Close the file
    MyFile.close();

}

string getFileName(int choice) //Test File Name module, 10 files in total present
{

	if(choice==1) return "eg1.txt";
	else if(choice==2) return "eg2.txt";
	else if(choice==3) return "eg3.txt";
	else if(choice==4) return "eg4.txt";
	else if(choice==5) return "eg5.txt";
	else if(choice==6) return "eg6.txt";
	else if(choice==7) return "eg7.txt";
	else if(choice==8) return "eg8.txt";
	else if(choice==9) return "eg9.txt";
	else if(choice==10) return "eg10.txt";
	else
	{
		cout<<"\n\nSpecified Test File NOT Present.\n\n";
		return "-";
	}
}
 
/* Function to print linked list */
void printLL(TopologicalOrder *head)
{
        TopologicalOrder* temp = head;

		cout<<"\n\n DFS Order: ";

        while (temp != NULL) 
		{
			cout << temp->vertex << " ";
            temp = temp->next;
        }
 }




// Utility function to print Adjacency List in file
void printInFile(DiaGraph *g,int V,string filename,int flag)
{
	filename.insert(0,"./Graphviz/");

	if(!flag)
	{
		ofstream MyFile(filename);

		// Write to the file
		MyFile << "digraph G{\n";

		for (int i = 0; i < V; i++)
		{
				// display adjacent vertices of vertex i
				adjNode* ptr=g->head[i];

				while (ptr != NULL) 
				{
					MyFile << i << "->" << ptr->val <<";\n";
					ptr = ptr->next;
				}

				if(ptr==NULL)
					MyFile << i <<";\n";

		}

		MyFile << "}";

		// Close the file
		MyFile.close();
	 }

	else
	{		
		ofstream MyFile(filename);

		// Write to the file
		MyFile << "digraph G{\n";

		for (int i = 0; i < V; i++)
		{
			// display adjacent vertices of vertex i
			adjNode* ptr=g->head[i];

			while (ptr != NULL) 
			{
				MyFile << i << "->" << ptr->val << " [style=bold,label="<<ptr->cost<<"]" <<";\n";
				ptr = ptr->next;
			}

			if(ptr==NULL)
				MyFile << i <<";\n";

		}

		MyFile << "}";

		// Close the file
		MyFile.close();
	}
}

// Utility function to print all adjacent vertices of given vertex
void display_AdjList(adjNode* ptr, int i,int flag)
{
	if(flag)
    {
		while (ptr != NULL) 
		{
			cout << "(" << i << ", " << ptr->val << ", " << ptr->cost << ") ";

			ptr = ptr->next;
		}

		cout << endl;
	}

	else
	{
		while (ptr != NULL) 
		{
			cout << "(" << i << ", " << ptr->val << ") ";

			ptr = ptr->next;
		}

		cout << endl;
	}
}



// Utility function for DFS Printing
void printDFS(int u,DiaGraph *g,int V,vector<int> pre,vector<int> post)
{
	string filename="DFS.gv";

	filename.insert(0,"./Graphviz/");

	ofstream MyFile(filename);

    // Write to the file
    MyFile << "digraph G{\n";

	vector<adjNode*> v; 

	for(int i=0;i<V;i++)
	{
		MyFile << i << "[style=bold,label=<<FONT COLOR=\"BLUE\" POINT-SIZE=\"24.0\">"<<i<<"</FONT><BR/>"<<pre[i]<<"/"<<post[i]<<">];\n";
	}
	
	for (int i = u; i < V; i++)
    {
        // display adjacent vertices of vertex i
        adjNode* ptr=g->head[i];

		 while (ptr != NULL) 
		 {
			 if(ptr->type=='t')
			 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"tree"<<"]" <<";\n";
				 

			 else if(pre[i]<pre[ptr->val] && post[i]>post[ptr->val])
			 {
				 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"forward"<<"]" <<";\n";
				 ptr->type='f';
			 }


			 else if(pre[i]>pre[ptr->val] && post[i]<post[ptr->val])
			 {
				 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"back"<<"]" <<";\n";
				 ptr->type='b';
			 }

			 else if(pre[i]>pre[ptr->val] && post[i]>post[ptr->val])
			 {
				 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"cross"<<"]" <<";\n";
				 ptr->type='c';
			 }

			 ptr = ptr->next;
		 }

    }



	for (int i = u-1; i >= 0; i--)
    {
        // display adjacent vertices of vertex i
        adjNode* ptr=g->head[i];

		 while (ptr != NULL) 
		 {
			 if(ptr->type=='t')
			 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"tree"<<"]" <<";\n";


			 else if(pre[i]<pre[ptr->val] && post[i]>post[ptr->val])
			 {
				 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"forward"<<"]" <<";\n";
				 ptr->type='f';
			 }


			 else if(pre[i]>pre[ptr->val] && post[i]<post[ptr->val])
			 {
				 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"back"<<"]" <<";\n";
				 ptr->type='b';
			 }

			 else if(pre[i]>pre[ptr->val] && post[i]>post[ptr->val])
			 {
				 MyFile << i << "->" << ptr->val << " [style=bold,label="<<"cross"<<"]" <<";\n";
				 ptr->type='c';
			 }

			 ptr = ptr->next;
		 }

    }



	// Number of nodes in graph
    /*for(int i = 0; i < V; i++)
        cout << "Node " << i << " Pre number "<< pre[i] << " Post number "<< post[i] << endl;*/

	MyFile << "}";

	// Close the file
    MyFile.close();
}

// Utility Function to check if element present in Vector or not
bool InVector(int num, vector<int> arr)
{
	for(int i=0;i<arr.size();i++)
			if(arr[i]==num)
				return true;

	return false;
}

// Utility Function to print graphs according to Tarjans Algoritm
void TarjanPrint(vector<vector<int>> scc,DiaGraph *g, int V,int flag)
{
	
	if(flag)
	{	
		
		ofstream MyFile("./Graphviz/TarjanSCC.gv");

		// Write to the file
		MyFile << "digraph G{\n";

		for (int i = 0; i < V; i++)
		{
			// display adjacent vertices of vertex i
			adjNode* ptr=g->head[i];

			while (ptr != NULL) 
			{
				for(int j=0;j<scc.size();j++)
				 {
					 if(scc[j].size()==1)
						  MyFile << scc[j][0] << " [style=bold]" <<";\n";

					 else if(InVector(i,scc[j]) && InVector(ptr->val,scc[j]))
					 {
						 MyFile << i << "->" << ptr->val << ";\n";
					 }
				}

				 ptr = ptr->next;
			}

		}

		MyFile << "}";

		// Close the file
		MyFile.close();
	}
}

// Utility Function to check edge between 2 components.
bool edgeBetween(vector<int> a,vector<int> b,vector<graphEdge> edges)
{
	for(int k=0;k<a.size();k++)
		{
			for(int l=0;l<b.size();l++)
			{
				for(int e=0;e<edges.size();e++)
				{
					if(edges[e].start_ver==a[k] && edges[e].end_ver==b[l])
					{
						return true;
					}
				}
				
			}
		}

	return false;
}

// Utility Function to check edge between 2 vertices
bool edgeBetween(int a,int b,vector<graphEdge> edges)
{
	for(int e=0;e<edges.size();e++)
	{
		if(edges[e].start_ver==a && edges[e].end_ver==b)
		{
						return true;
		}
	}

	return false;
}


// Utility Function to check edge between 2 components & return edgeset
vector<int> edgeBetweenValues(vector<int> a,vector<int> b,vector<graphEdge> edges)
{
	vector<int> bw(3);

	for(int k=0;k<a.size();k++)
		{
			for(int l=0;l<b.size();l++)
			{
				for(int e=0;e<edges.size();e++)
				{
					if(edges[e].start_ver==a[k] && edges[e].end_ver==b[l])
					{
						bw[0]=a[k];
						bw[1]=b[l];
						bw[2]=edges[e].weight;
						return bw;
					}
				}
				
			}
		}
		
		
	return bw;

}



// A utility function to create a new Min Heap Node
MinHeapNode* newMinHeapNode(int v,int dist)
{
    MinHeapNode* minHeapNode = new MinHeapNode();
    minHeapNode->v = v;
    minHeapNode->dist = dist;

    return minHeapNode;
}

// A utility function to create a Min Heap
MinHeap* BuildMinHeap(int capacity)
{
    MinHeap* minHeap =new MinHeap();

    minHeap->pos = new int[capacity * sizeof(int)];
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->arr =new MinHeapNode*[capacity *sizeof(MinHeapNode*)];

    return minHeap;
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b)
{
    MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}
 
// A standard function to heapify at given idx This function also updates position of nodes when they are swapped. Position is needed for decreaseKey()
void minHeapify(MinHeap* minHeap,int index)
{
    int smallest, left, right;
    smallest = index;
    left = 2 * index + 1;
    right = 2 * index + 2;
 
    if (left < minHeap->size && minHeap->arr[left]->dist < minHeap->arr[smallest]->dist )
      smallest = left;
 
    if (right < minHeap->size && minHeap->arr[right]->dist < minHeap->arr[smallest]->dist )
      smallest = right;
 
    if (smallest != index)
    {
        // The nodes to be swapped in min heap
        MinHeapNode *smallestNode = minHeap->arr[smallest];
        MinHeapNode *idxNode = minHeap->arr[index];
 
        // Swap positions
        minHeap->pos[smallestNode->v] = index;
        minHeap->pos[idxNode->v] = smallest;
 
        // Swap nodes
        swapMinHeapNode(&minHeap->arr[smallest],&minHeap->arr[index]);
 
        minHeapify(minHeap, smallest);
    }
}
 
// A utility function to check if the given minHeap is empty or not
int isEmpty(MinHeap* minHeap)
{
    return minHeap->size == 0;
}
 
// Standard function to extract minimum node from heap
MinHeapNode* extractMin(MinHeap* minHeap)
{
    if (isEmpty(minHeap))
        return NULL;
 
    // Store the root node
    MinHeapNode* root = minHeap->arr[0];
 
    // Replace root node with last node
    MinHeapNode* lastNode =minHeap->arr[minHeap->size - 1];
    minHeap->arr[0] = lastNode;
 
    // Update position of last node
    minHeap->pos[root->v] = minHeap->size-1;
    minHeap->pos[lastNode->v] = 0;
 
    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);
 
    return root;
}
 
// Function to decrease dist value of a given vertex v. This function uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(MinHeap* minHeap, int v, int dist)
{
    // Get the index of v in  heap array
    int i = minHeap->pos[v];
 
    // Get the node and update its dist value
    minHeap->arr[i]->dist = dist;
 
    // Travel up while the complete tree is not hepified. This is a O(Logn) loop
    while (i && minHeap->arr[i]->dist <minHeap->arr[(i - 1) / 2]->dist)
    {
        // Swap this node with its parent
        minHeap->pos[minHeap->arr[i]->v] =(i-1)/2;
        minHeap->pos[minHeap->arr[(i-1)/2]->v] = i;
        swapMinHeapNode(&minHeap->arr[i], &minHeap->arr[(i - 1) / 2]);
 
        // move to parent index
        i = (i - 1) / 2;
    }
}
 
// A utility function to check if a given vertex 'v' is in min heap or not
bool isInMinHeap(MinHeap *minHeap, int v)
{
   if (minHeap->pos[v] < minHeap->size)
     return true;

   return false;
}
 
// A utility function used to print the solution of Dijikstra in File & console
void printArr(int dist[], int n,int par[],vector<graphEdge> edges,int t)
{
    cout<<"\n\nVertex   Distance from Source   Parent\n";
    for (int i = 0; i < n; ++i)
	{
		if(dist[i]==INT_MAX)
		cout<<i<<" \t\t Unreachable     "<<par[i]<<"\n";

		else if(par[i]==-1)
			cout<<i<<" \t\t "<<dist[i]<<"\t\t"<<"--"<<"\n";

		else
			cout<<i<<" \t\t "<<dist[i]<<"\t\t"<<par[i]<<"\n";


	}

	cout<<"\n\nShortest Path cost to "<<t<<" from given vertex is: ";
    for (int i = 0; i < n; ++i)
	{
		if(i==t)
		{
			if(dist[i]==INT_MAX)
				cout<<"Unreachable\n\n";

			else
				cout<<dist[i]<<"\n\n";
		}
	}

	ofstream MyFile("./Graphviz/Dijikstra.gv");

    // Write to the file
    MyFile << "digraph G{\n";

	for (int i = 0; i < n; ++i)
	{
		if(par[i]==-1)
			continue;

		else
		{
			int linkval;

			for(int e=0;e<edges.size();e++)
			{
				if(par[i]==edges[e].start_ver && i==edges[e].end_ver) // Print cost of parent------>child link
				{
					linkval=edges[e].weight;
					MyFile << par[i] << "->" << i << " [style=bold,label="<<linkval<<"]" <<";\n";
				}
			}
			
		}
	}
	
	MyFile << "}";

	// Close the file
    MyFile.close();

}

/********** Utility Functions End ********** */

/* *********Assignment Questons************* */


// Question-1 Standard DFS
void DFS(int u,DiaGraph *g, vector<int> &pre,vector<int> &post,vector<int> &vis, int V, int *Time,TopologicalOrder **list,int flag)
{
	
	// Storing the pre number whenever
    // the node comes into recursion stack
    pre[u] = *Time;

    // Increment time
    *Time=*Time+1;
    vis[u] = 1;

	if(flag)
	cout<<u<<" ";

	adjNode* ptr=g->head[u];

	while(ptr)
    {
		int v=ptr->val;

		if (vis[v] == 0)
		{
			ptr->type='t';
			DFS(v, g, pre, post, vis,V,Time,list,flag);
		}

		ptr=ptr->next;
	}
 
    // Storing the post number whenever
    // the node goes out of recursion stack
    post[u] = *Time;
	
	if(!(*list))
	{
		*list=new TopologicalOrder();
		(*list)->vertex=u;
		(*list)->next=NULL;
	}

	else
	{
		TopologicalOrder* temp=new TopologicalOrder();

		temp->vertex=u;
		temp->next=*list;

		*list=temp;
	}

    *Time=*Time+1;

}

// Question-1 Standard DFS which also helps to find Topological order
TopologicalOrder* DFS(int u,DiaGraph *g, vector<int> &pre,vector<int> &post,vector<int> &vis, int V,int flag)
{
	TopologicalOrder* list=NULL;

	int Time = 1;
	DFS(u,g,pre,post,vis,V,&Time,&list,flag);

	for(int i=0;i<V;i++)
	{
		if(vis[i]==0)
			DFS(i, g, pre, post, vis,V,&Time,&list,flag);
	}

	return list;
}

// Utility Function Return Topological Ordering 
TopologicalOrder* Topological_DAG(DiaGraph *g,int V)
{
	TopologicalOrder* list=NULL;

	int flag=0;

	vector<int> pre(V);
    vector<int> post(V);
 
    // Visited array
    vector<int> vis(V);

	list=DFS(0, g, pre, post, vis,V,flag);

	return list;
}


// Utility Function for DFS Traversal in Tarjans Algoritm
void DFS_Tarjan(int u,vector<int> &df,vector<int> &low,stack<int> &stk, vector<bool> &InStk,DiaGraph *g,int V,int *t,vector<vector<int>> &scc,int flag)
{
	
	df[u]=low[u]=*t;

	*t=*t+1;

	stk.push(u);

	InStk[u]=true;

	adjNode* ptr=g->head[u];

	while(ptr)
    {
		int v=ptr->val;

		if (df[v] == 0) 
		{
			
			DFS_Tarjan(v, df, low, stk, InStk,g,V,t,scc,flag);
			low[u]=min(low[u],low[v]); // tree edge
		}

		else if(InStk[v]) // back edge
		{
			low[u]=min(low[u],df[v]);
		}

		ptr=ptr->next;
	}
 
    if(low[u]==df[u])
	{
		if(flag)
		cout<<"\n\nSCC: ";

		vector<int> temp;

		while(stk.top() != u)
		{
			if(flag)
			cout<<stk.top()<<" ";
			temp.push_back(stk.top());
			InStk[stk.top()]=false;
			stk.pop();
		}

		if(flag)
		cout<<stk.top()<<"\n";
		temp.push_back(stk.top());
		InStk[stk.top()]=false;
		stk.pop();

		
		scc.push_back(temp);
	}

}

// Utility Function for DFS Traversal in Tarjans Algoritm
vector<vector<int>> DFS_Tarjan(int u,vector<int> &df,vector<int> &low,stack<int> &stk, vector<bool> &InStk,DiaGraph *g,int V,int flag)
{
	int t=1;

	vector<vector<int>> scc;

	DFS_Tarjan( u,df,low,stk, InStk, g, V,&t,scc,flag);

	

	return scc;
}


// Question-2 Solution
vector<vector<int>> Tarjans_SCC(DiaGraph *g,int V, int flag)
{
	vector<int> df(V,0);
	vector<int> low(V,0);
	vector<vector<int>> scc;
	vector<vector<int>> tmp;

	vector<bool> InStk(V,false);

	stack<int> stk;

	int i=0;

	for(;i<V;i++)
	{
		if(df[i]==0)
		{
			tmp=scc;

			scc=DFS_Tarjan(i,df,low,stk,InStk,g,V,flag);

			if(scc.size()==1)
			{
				vector<int> v=scc[scc.size()-1];
				scc=tmp;
				scc.push_back(v);
			}

		}
	}
	
	TarjanPrint(scc,g,V,flag);

	return scc;
}


// Question-5 Solution
void Dijikstra(DiaGraph* graph,int s,int t,int V,vector<graphEdge> edges)
{

	// dist values used to pick minimum weight edge 
    int *dist=new int[V];    
	int *parent=new int[V];
 
    // minHeap represents set E
    MinHeap* minHeap = BuildMinHeap(V);
 
    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v)
    {
        dist[v] = INT_MAX;
		parent[v]=-1;
        minHeap->arr[v] = newMinHeapNode(v,dist[v]);
        minHeap->pos[v] = v;
    }
 
    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->arr[s] = newMinHeapNode(s, dist[s]);
    minHeap->pos[s]   = s;
    dist[s] = 0;
	parent[s]=-1;
    decreaseKey(minHeap, s, dist[s]);
 
    // Initially size of min heap is equal to V
    minHeap->size = V;
 
    // In the following loop, min heap contains all nodes who are NOT closed
    while (!isEmpty(minHeap))
    {
        // Extract the vertex with minimum distance value
        MinHeapNode* minHeapNode = extractMin(minHeap);
       
        // Store the extracted vertex number
        int u = minHeapNode->v;
 
        // Traverse through all adjacent vertices of u (the extracted vertex) and update their distance values
        adjNode* ptr = graph->head[u];

        while (ptr != NULL)
        {
			int v = ptr->val;
 
            // If shortest distance to v is not finalized yet, and distance to v through u is less than its previously calculated distance
            if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX && ptr->cost + dist[u] < dist[v])
            {
                dist[v] = dist[u] + ptr->cost;

				parent[v]=u;
 
                // update distance value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }
            ptr = ptr->next;
        }
    }
 
    // print the calculated shortest distances
    printArr(dist, V,parent,edges,t);

}




// Question-4 Solution
bool IsSemiConnected(DiaGraph *diagraph,int V,vector<graphEdge> edges)
{
	TopologicalOrder* list=NULL;

	int flag=0;

	// Step-1 Find SCC

	vector<vector<int>> scc=Tarjans_SCC(diagraph,V,flag);
	
	int V2=scc.size(); // Number of Vertices in new Component graph 

	int E2=0;

	vector<graphEdge> newedges;

	bool edgeExist;

	if(V2==1) // Only 1 Component means Strongly Connected graph
		return true;


	// Step-2 Find Component Graph

	for(int i=0; i <scc.size(); i++)
	{
		for(int j=i+1; j <scc.size(); j++)
		{
			if(i != j)
			{
				edgeExist=edgeBetween(scc[i],scc[j],edges);

				if(edgeExist)
				{
					graphEdge g;

					g.start_ver=i;
					g.end_ver=j;
			
					newedges.push_back(g);
					E2++;
				}


				edgeExist=edgeBetween(scc[j],scc[i],edges);

				if(edgeExist)
				{
					graphEdge g;

					g.start_ver=j;
					g.end_ver=i;
			
					newedges.push_back(g);
					E2++;
				}
			}
		}
	}

	DiaGraph *graph2=new DiaGraph(newedges, E2, V2);

	/*// print adjacency list representation of graph
    cout<<"Graph adjacency list "<<endl<<"(start_vertex, end_vertex, weight):"<<endl;

    for (int i = 0; i < V2; i++)
    {
        // display adjacent vertices of vertex i
        display_AdjList(graph2->head[i], i);
    }
	*/

	// Step-3 Find Topological Ordering of Component Graph
	list=Topological_DAG(graph2,V2);

	TopologicalOrder *ptr;

	ptr=list;

	int count=0;


	// Step-4 For Topological Ordering check if there is Edge between Vi and Vi+1
	while(ptr) 
	{
		if(ptr->vertex==0)
			count++;

		else
		{
			for(int e=0;e<newedges.size();e++)
			{
				if(newedges[e].start_ver==ptr->vertex && newedges[e].end_ver==ptr->vertex-1)
				{
					count++;
					break;
				}
			}
		}

		ptr=ptr->next;
	}

	// If there is Edge between every Vi and Vi+1 then return true
	if(count == V2)
		return true;

	else
	return false;

}


// Question-3 Solution:-
void minEdges(DiaGraph *diagraph,int V,vector<graphEdge>edges)
{

	vector<vector<int>> scc=Tarjans_SCC(diagraph,V,0);
	
	int V2=V; // Number of Vertices in new graph must be same as per requirement

	int E2=0;

	vector<graphEdge> newedges;

	bool edgeExist;

	for(int i=0; i <scc.size(); i++)
	{

		if(scc[i].size()>1)
		{
			for(int k=0;k<scc[i].size()-1;k++)
			{
				for(int m=k+1;m<scc[i].size();m++)
				{
					edgeExist=edgeBetween(scc[i][k],scc[i][m],edges);

					if(edgeExist)
					{
						graphEdge g;

						g.start_ver=scc[i][k];
						g.end_ver=scc[i][m];
			
						newedges.push_back(g);
						E2++;
					}

					edgeExist=edgeBetween(scc[i][m],scc[i][k],edges);

					if(edgeExist)
					{
						graphEdge g;

						g.start_ver=scc[i][m];
						g.end_ver=scc[i][k];
			
						newedges.push_back(g);
						E2++;
					}
				}
			}
		}

		for(int j=i+1;j<scc.size();j++)
		{
			edgeExist=edgeBetween(scc[i],scc[j],edges);

			if(edgeExist)
			{
				vector<int> bw(3);
				bw=edgeBetweenValues(scc[i],scc[j],edges);

				graphEdge g;

				g.start_ver=bw[0];
				g.end_ver=bw[1];
				g.weight=bw[2];
			
				newedges.push_back(g);
				E2++;
			}

			edgeExist=edgeBetween(scc[j],scc[i],edges);

			if(edgeExist)
			{
				vector<int> bw(3);
				bw=edgeBetweenValues(scc[j],scc[i],edges);

				graphEdge g;

				g.start_ver=bw[0];
				g.end_ver=bw[1];
				g.weight=bw[2];
			
				newedges.push_back(g);
				E2++;
			}

		}
	}

	DiaGraph *graph=new DiaGraph(newedges, E2, V2);

	printInFile(graph,V,"MinEdge_Graph.gv",0);

	// print adjacency list representation of graph
    cout<<"\n\nMinEdge Graph adjacency list of form (start_vertex, end_vertex):"<<endl;

    for (int i = 0; i < V; i++)
    {
        // display adjacent vertices of vertex i
        display_AdjList(graph->head[i], i,0);
    }


	return;
}

/* ***** Main Function ****** */

int main()
{
	char c;
	int choice,num;

	cout<<"\n\n\n********* Menu Driven Mode *********\n\n\n";

	while(1)
	{
		cout<<"\n\nDo You want to create your own testfile or use pre-build testfile which is already present in input folder ?\n\n";
		cout<<"\n\nWrite 1 to create your own file or 0 for using already made ten test files: ";
		cin>>choice;

		string input,filename;

		if(choice == 1)
		{
			cout<<"\n\nEnter any name of file you want [For example: test.txt,check.txt etc.]          (Please ensure to write .txt at the end) : ";
			cin>>input;

			int len =input.length();

			if(input[len-1]!='t' ||input[len-2]!='x' ||input[len-3]!='t' || input[len-4]!='.')
			{
				cout<<"\n\nYou entered filename without .txt at the end !!!\n\n";
				cout<<"Try again [y/n] ?: ";
				cin>>c;

				if(c=='y')
				continue;

				else
					continue;
			}

			createOwnFile(input);
			filename=input;

			
		}


		else if(choice == 0)
		{
			cout<<"\n\nEnter Test_File number [ Between 1-10 only]: ";
			cin>>choice;

			filename=getFileName(choice);

			if(filename=="-")
				continue;
			
		}

		else
		{
			cout<<"\n\nInvalid Choice!!\n\n";
			continue;
		}

		//char *filename="eg2.txt";

		if(filename!="-")
		{
			int V,E; // Vertex start from 0,1,2......

			vector<graphEdge> edges;

			string myText;

			filename.insert(0,"./input/");

			ifstream MyReadFile(filename);

			int line=1;

			while (getline (MyReadFile, myText))
			{
				int len=myText.length();

				int i=0;

				if(len<3) // If no text in line
				break;

				if(line==1) // line 1 must have 2 numbers only
				{
					string str1;
					string str2;

					while(myText[i] != ' ')
					{
						str1+=myText[i];
						i++;
					}

					stringstream tonum1(str1);

					tonum1 >> V;

					i++;

					while(i < len)
					{
						str2+=myText[i];
						i++;
					}

					stringstream tonum2(str2);

					tonum2 >> E;
					line++;
				}

				else // Rest line must have 3 numbers only
				{
					string str1;
					string str2;
					string str3;

					int num1,num2,num3;

					while(myText[i] != ' ')
					{
						str1+=myText[i];
						i++;
					}

					stringstream tonum1(str1);

					tonum1 >> num1;

					i++;

					while(myText[i] != ' ')
					{
						str2+=myText[i];
						i++;
					}

					stringstream tonum2(str2);

					tonum2 >> num2;

					i++;

					int flag=0;

					if(myText[i]=='-') // Negative Weight
					{
						flag=1;
						i++;
					}	

					while(i != len)
					{
						str3+=myText[i];
						i++;
					}

					stringstream tonum3(str3);

					tonum3 >> num3;

					if(flag==1)
						num3=-(num3);

					addEdge(edges,num1,num2,num3); // Add Edge to edge set

					line++;
				}
			}

			DiaGraph *diagraph=new DiaGraph(edges, E, V); // Use edge set to initialise Directed graph

			// print adjacency list representation of graph
			cout<<"\n\nAdjacency list of Graph File of form "<<"(start_vertex, end_vertex, weight) is:\n"<<endl;

			for (int i = 0; i < V; i++)
			{
				// display adjacent vertices of vertex i
				display_AdjList(diagraph->head[i], i,1);
			}
	

			printInFile(diagraph,V,"originalGraph.gv",1);

			cout<<"\n\nGraphviz File(DOT file) is stored in the Graphviz Directory as .gv file with name: originalGraph.gv. Please Check & run it as per readme file.\n\n";

			bool endLoop=false;

			while(!endLoop)
			{
				cout<<"\n\nEnter Number corresponding to operation-\n\n1. DFS Traversal on Graph\n\n2. SCC Using Tarjan's algorithm\n\n3. Minimum Edge Graph\n\n4. Check If Graph Semi-connected\n\n5. Dijikstra Shortest Path\n\n6. Exit\n\n";

				cout<<"Choice: ";

				cin>>choice;

				switch(choice)
				{
					case 1:
					{
						/* ****** DFS Start ****** */

						int start_dfs;

						cout<<"\n\nEnter Start vertex of DFS traversal [between "<<0<<" - "<<V-1<<"]: ";
						cin>>start_dfs;
						
						if(start_dfs>V-1)
						{
							cout<<"Vertex out of limits Specified";
							continue;
						}

						cout<<"\n\nDFS Order: ";

						vector<int> pre(V);
						vector<int> post(V);
 
						// Visited array
						vector<int> vis(V);

						TopologicalOrder* o=DFS(start_dfs, diagraph, pre, post, vis,V,1);

						printDFS(start_dfs,diagraph,V,pre,post);

						
						cout<<"\n\n";
						
						/* **** DFS End **** */

						cout<<"\n\nGraphviz File(DOT file) is stored in the Graphviz Directory as .gv file with name: DFS.gv. Please Check & run it as per readme file.\n\n";

						break;
					}

					case 2:
					{
						Tarjans_SCC(diagraph,V,1);

						cout<<"\n\nGraphviz File(DOT file) is stored in the Graphviz Directory as .gv file with name: TarjanSCC.gv. Please Check & run it as per readme file.\n\n";

						break;
					}

					case 3:
					{
						/* MinEdges */
						minEdges(diagraph,V,edges);

						cout<<"\n\nGraphviz File(DOT file) is stored in the Graphviz Directory as .gv file with name: MinEdge_Graph.gv. Please Check & run it as per readme file.\n\n";

						break;
					}

					case 4:
					{
						/*Semiconnected*/
						bool res=IsSemiConnected(diagraph,V,edges);

						if(res)
						cout<<"\n\nGiven Graph is Semi-connected\n\n";

						else
						cout<<"\n\nGiven Graph is NOT Semi-connected\n\n";

						

						break;
					}

					case 5:
					{
						/* Dijikstra Shortest Path*/
						int s;
						int t;

						cout<<"\n\nSource[Please enter no. between "<<0<<"-"<<V-1<<" only]: ";
						cin>>s;

						cout<<endl;

						cout<<"\n\nDestination[Please enter no. between "<<0<<"-"<<V-1<<" only]: ";
						cin>>t;

						if(s <=-1 || s > V-1 || t <=-1 || t > V-1)
						{
							cout<<"\n\nYou entered something wrong. Do this step Again\n\n";
							break;
						}

						Dijikstra(diagraph,s,t,V,edges);

						cout<<"\n\nGraphviz File(DOT file) is stored in the Graphviz Directory as .gv file with name: Dijikstra.gv. Please Check & run it as per readme file.\n\n";

						break;
					}

					case 6:
						endLoop=true;
						break;


					default:
						printf("\nInvalid choice\n");

				}
			}

	

		}

		cout<<"\n\nDo you want to Test other test file?[y/n] : ";
		cin>> c;

		if(c=='y')
			continue;

		else
		{
			cout<<"\n\n*******Menu Driven Mode exiting*******\n\n";
			exit(0);
		}

	}




	return 0;
}

