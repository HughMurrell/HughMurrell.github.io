# usage
# julia urls.txt > clean_urls.txt
# xargs -n 1 curl -O < clean_urls.txt
#

open("urls.txt") do f
 
  # read till end of file
  while ! eof(f)  
 
     # read a new / next line for every iteration           
     s = readline(f)          
     url=split(s,'"');
     if length(url) > 1	 
     	println("https://www.pgnmentor.com/$(url[2])")
     end
  end
 
end
