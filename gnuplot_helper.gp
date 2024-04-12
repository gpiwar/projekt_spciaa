set grid 
set terminal gif animate delay 0.1
set output output_gif

index=0
while (index<=last_index){
	fname = sprintf('%s_%d.data', input_prefix, index)
	plot fname with image
	index=index+step
} 
