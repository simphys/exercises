reset
set grid
set datafile commentschar "#@"
unset key
set term pdf enhanced
set xlabel "Phi"
set ylabel "Psi"
set output "./report/plt/rama.pdf"
plot "./alanine_dipeptide/rama_alaninedipeptide.xvg" u 1:2 w p pt 13
unset term
