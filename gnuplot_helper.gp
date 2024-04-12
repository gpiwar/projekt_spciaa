unset border
unset xtics
unset ytics
unset raxis
unset rtics
unset key
unset colorbox

set lmargin at screen 0
set rmargin at screen 1
set tmargin at screen 0
set bmargin at screen 1

set terminal gif animate delay 0.1
set output output_gif

if (input_prefix eq 'out'){
	set cbrange [0:1400]
}

index=0
while (index<=last_index){
	fname = sprintf('%s_%d.data', input_prefix, index)
	plot fname with image title sprintf('%s%d.data', input_prefix, index)
	index=index+step
}
