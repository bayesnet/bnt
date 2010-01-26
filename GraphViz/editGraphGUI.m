function g
%here is how one creates a function ("callback") which does something
%(prints the node label) when you click on the node's text in Matlab figure.
%
% Leon Peshkin  http://www.ai.mit.edu/~pesha
%
%draw_graph(...)

     % "gca" is the current "axes" object, parent of all objects in figure
     % "gcbo" is the handle of the object whose callback is being executed
     % "findall" gives handles to all elements of a given type in the figure
text_elms = findall(gca,'Type','text');  
for ndx = 1:length(text_elms)
  callbk = 'my_call(str2num(get(gcbo,''String'')))'; 
  set(text_elms(ndx), 'ButtonDownFcn', callbk);  % assume the node label is a number
end
