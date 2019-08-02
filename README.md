# New_only_forLc_Dfinder
1) First, delet all the useless things:  For Lc, we do not need "Dfinder::TkCombinationResFast" (https://github.com/taweiXcms/Bfinder/blob/Dfinder/Bfinder/src/Dfinder.cc), and we do not need resonance particles, so we delet all of these things. This will help you speed up a little.

2) The first biggest change is "TkCombinationPermutation_Lc_v3" (this is changed from "TkCombinationPermutation" in the above github link.) In this part, we change "mass cut" to "mass square cut" (line 1622). Also we change the pT cut, rapidity cut (line 1617, 1624.) In additon doing some transformations of the cut variables, we also defind "TrackXYZP2" before this part. Here we record the charge of tracks to match the permutation of reconstruction Lc. We use 4(q1+1)+2(q2+1)+(q3+1). For each specific particle, this needs to be changed. [NOTE: THIS PART IS ONLY WORKS FOR Lc]

3) This second biggest change is from line 1741-1758, here have to match the permutation with the track mass.

4) Since we change the way to find the permutation (step2 and step 3), We have to change this part: line 985, we do not need to loop all the permutations, so we comment line 985. 

5) Combine step 2,3 and 4 will help speeding up a lot. (The reason is that: we loop all the permutaions previously, and now we go through once and find this is whether Lc permutation or not.)
