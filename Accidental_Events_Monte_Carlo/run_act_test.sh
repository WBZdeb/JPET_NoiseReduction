#!/bin/bash

output_file="act_output.txt"

B_values=(10000.0 30000.0 60000.0 100000.0 300000.0 600000.0 1000000.0 3000000.0 6000000.0 8000000.0 9000000.0 10000000.0 30000000.0 20000000.0 40000000.0 50000000.0)
W_values=(2000 1000 1000 500 250 100 50 20 10 5 5 2 2 2 2 2)

pattern="All randoms:"

for i in "${!B_values[@]}"; do
	B="${B_values[i]}"
	W="${W_values[i]}"

	echo "Iteration $((i + 1)): Running with B = $B and W = $W"

	echo -n "$B " >> "$output_file"
	
	root -l -q "testDTW.cpp($W, $B)" 2>&1 | grep "$pattern" | awk '{print $3}' >> "$output_file"
	
	#echo -e "\n" >> "$output_file"
done

echo "All done."
