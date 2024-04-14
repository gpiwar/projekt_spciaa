set grid 
set terminal gif animate delay 0.1
set output 'oil_out.gif'
set cbrange[0:0.0012]

n=0
while (n<50000){
	fname = sprintf('out_%d.data', n)
	plot fname with image
	n=n+100
} 

set output
