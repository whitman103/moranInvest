set terminal gif animate delay 5
set output "cancerGrowth.gif"

set xrange[15:35]
set yrange[15:35]
set zrange[0:15]

do for [i=0:199]{
	set view 49, i
	splot "dataBin//cancerCells_".i.".txt" u 1:2:3:4 with points ps 1.5 pt 7 lc variable
}

unset output