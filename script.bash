cmake . && make
./gen 1 0 1
./solution preprocess --nodes n_gen.txt --edges e_gen.txt --output g.txt
./solution search --graph g.txt --input q_gen.txt --output res.txt --full-output
./validator n_gen.txt e_gen.txt q_gen.txt res.txt
