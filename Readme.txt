***How To Run .cpp FILE***

Pre-Requisites: 

Please ensure that "214101020_Graphs.cpp", "input folder" and "Graphviz folder" are present within SAME directory. 
NOTE : Placing them in same place is IMPORTANT, otherwise program will give Segment Fault.
         
For Example:
			1. Create a "CPP" named folder in easily accessible location. In my machine, I placed it in D:/ drive.
			2. Inside CPP folder, place "214101020_Graphs.cpp", "input folder" & "Graphviz folder".

Running commands:

1. On Linux Machine use-
			(i) Use following command after traversing to CPP folder (implicitly in cmd/shell).

				CPP/> g++ 214101020_Graphs.cpp
				CPP/> ./a.out

			

2. On Windows you may use IDE to directly run or use Cygwin to simulate Linux like environment, And then use above commands.

3. You will be shown following outputs-

		
				*** MENU DRIVEN MODE ***

				Do You want to create your own testfile or use pre-build testfile which is already present in input folder ?
  				
				Write 1 to create your own file or 0 for using already made ten test files:
				
				

				*If you press 0, there are 10 Distinct Test_files namely eg1,eg2,eg3.....eg10 in input folder. Choose any 1 by entering number between 1-10.

				You will be shown an Adjacency list and 5 Assignment Solutions.

				Evaluate Test file graph's DFS, Strongly Connected Components, Shortest paths etc. by pressing suitable options.






				***NOTE: If you press 1, Enter the data very carefully , otherwise error may turn up.
				
				After pressing 1, First you need to specify filename with .txt extension. Enter any name of file you want [For example: test.txt,check.txt etc.].
				
				***Having .txt as extension is IMPORTANT. So ensure that you write it at end of file name [For example: test.txt,check.txt etc. are valid filenames].

				After this step, you need to specify details like no. of vertices=V, edges=E. 

				***Based on no. of vertices, 0,1,2,3......V-1 named vertices will be created. Based on Edges E, you need to specify a directed edge from u---->v with weight w.

				On Successfully creating a file in "input directory". You will be shown an Adjacency list and 5 Assignment Solution.

				Evaluate your graph's DFS, Strongly Connected Components, Shortest paths etc. by pressing suitable options.



				

				The Ouput of each of the 5 Operations can be seen visually & on Console too. 

				Five Output files will be present in "Graphviz folder". Run each of them as follows to see Visual Graph.
				
				
				
				
				
			
4. IMPORTANT NOTE: For visualising the .gv file of "214101020_Graphs.cpp" program, Please take absolute precautions. Note that after each operation, a DOT file will be generated in Graphviz folder.
    


Steps to visualize "originalGraph.gv", "DFS.gv", "TarjanSCC.gv", "MinEdge_Graph.gv" & "Dijikstra.gv" in Windows:

     A) Ensure that graphviz is installed in your machine.

     B) Go to cmd. And then to CPP folder location. In my machine I placed CPP folder in D drive. Then to Graphviz folder
	 
                                                          ---------->    D:\CPP\        
														  ---------->    D:\CPP\Graphviz\
                                                         

      C) Write command in following format-

                                    D:\CPP\Graphviz\> "dot executable location in your computer" -Tpng <filename.gv> -o <newfilename.png>

   		        For above example in my pc :

                            						-----> D:\CPP\Graphviz\> "C:\Program Files\Graphviz2.38\bin\dot.exe" -Tpng originalGraph.gv -o originalGraph.png
							-----> D:\CPP\Graphviz\> "C:\Program Files\Graphviz2.38\bin\dot.exe" -Tpng DFS.gv -o DFS.png 
							-----> D:\CPP\Graphviz\> "C:\Program Files\Graphviz2.38\bin\dot.exe" -Tpng TarjanSCC.gv -o TarjanSCC.png
							-----> D:\CPP\Graphviz\> "C:\Program Files\Graphviz2.38\bin\dot.exe" -Tpng MinEdge_Graph.gv -o MinEdge_Graph.png 
							-----> D:\CPP\Graphviz\> "C:\Program Files\Graphviz2.38\bin\dot.exe" -Tpng Dijikstra.gv -o Dijikstra.png 

D) Open png image to see output.



Steps to visualize "originalGraph.gv", "DFS.gv", "TarjanSCC.gv", "MinEdge_Graph.gv" & "Dijikstra.gv" in Linux:

	A) Ensure that graphviz is installed in your machine.
	
	B) Similar command works in linux too but instead of "C:\Program Files\Graphviz2.38\bin\dot.exe" only mention dot. Example
	
							-----> CPP\Graphviz\> dot -Tpng originalGraph.gv -o originalGraph.png
							-----> CPP\Graphviz\> dot -Tpng DFS.gv -o DFS.png 
							-----> CPP\Graphviz\> dot -Tpng TarjanSCC.gv -o TarjanSCC.png
							-----> CPP\Graphviz\> dot -Tpng MinEdge_Graph.gv -o MinEdge_Graph.png 
							-----> CPP\Graphviz\> dot -Tpng Dijikstra.gv -o Dijikstra.png 
							
	C) Open png image to see output.
