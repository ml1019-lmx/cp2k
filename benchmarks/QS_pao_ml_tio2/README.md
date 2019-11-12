# Quickstep Polarized Atomic Orbital - TiO2

Large scale benchmark for PAO ML (and/or LS DFT in general).

## How to Run the Benchmark

Bunzip2 all files, have [`BASIS_MOLOPT`](../../data/BASIS_MOLOPT) and [`GTH_POTENTIALS`](../../data/GTH_POTENTIALS) available (from [cp2k/data](../../data/).

For tuning purposes, the length of the full benchmark can be reduced in the following ways:

- by reducing the number of MD steps (STEPS 20 -> STEPS 5)
- by doing only an energy calculation (RUN_TYPE MD -> RUN_TYPE ENERGY)
- by doing an energy calculation with few SCF steps (MAX_SCF 50 -> MAX_SCF 5)

## Results Archive

Reference energies and timings using CP2K svn:17405, Piz Daint, Cray XC30, 1024 nodes, SB + K20X, +- 8Gb mem used per node

### Output File

(last column timings per step)

```
> grep SCF out | head -n 10
 ------------------------------ Linear scaling SCF -----------------------------
 SCF     1   -552962.443966503   -552962.443966503   28.188702
 SCF     2   -553391.806569889      -429.362603386   32.410283
 SCF     3   -553756.740393045      -794.296426542   12.872345
 SCF     4   -554517.794343487     -1555.350376984   13.497826
 SCF     5   -554728.247042532      -210.452699045   33.469438
 SCF     6   -554849.946612308      -332.152268821   13.510741
 SCF     7   -554873.984235220      -356.189891734   17.914532
 SCF     8   -554915.784987915       -41.800752695   38.287248
 SCF     9   -554918.728930513       -44.744695293   14.047717
```

### MD Energy File

```
> cat pao_ml_md-1.ener
#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]
         0            0.000000       110.495412573       300.000000000   -554937.745690580   -554827.250278007         0.000000000
         1            0.500000       106.078997716       288.009235623   -554933.183067304   -554827.104069588      1606.181380667
         2            1.000000       102.069290733       277.122701359   -554929.060740767   -554826.991450034       437.939860202
         3            1.500000       100.198165911       272.042513562   -554927.153307400   -554826.955141488       449.515135792
         4            2.000000       101.022853760       274.281578054   -554928.022891735   -554827.000037975       375.188747253
         5            2.500000       103.892195371       282.071969193   -554930.981726713   -554827.089531342       378.129010544
         6            3.000000       107.422614470       291.657215360   -554934.608040349   -554827.185425879       363.554193239
         7            3.500000       110.115891922       298.969584414   -554937.373128860   -554827.257236939       365.110804793
         8            4.000000       110.802924876       300.834909692   -554938.098364813   -554827.295439937       367.848488379
         9            4.500000       108.870951325       295.589514868   -554936.164765656   -554827.293814332       366.621533569
        10            5.000000       104.372470458       283.375937591   -554931.619703179   -554827.247232721       371.757427932
        11            5.500000        98.092255253       266.324871693   -554925.248720303   -554827.156465049       366.436341383
        12            6.000000        91.519001482       248.478192945   -554918.555918355   -554827.036916873       438.594749461
        13            6.500000        86.552351413       234.993515291   -554913.479003540   -554826.926652127       363.525124027
        14            7.000000        84.878892773       230.450000040   -554911.758150807   -554826.879258034       361.602635782
        15            7.500000        87.178039466       236.692286411   -554914.109975162   -554826.931935696       362.720366624
        16            8.000000        92.619988766       251.467422791   -554919.692436207   -554827.072447441       360.480368018
        17            8.500000        99.082357198       269.013042870   -554926.321931528   -554827.239574330       362.877061196
        18            9.000000       104.007271910       282.384407159   -554931.375035372   -554827.367763462       364.628673398
        19            9.500000       105.405738067       286.181305486   -554932.815285468   -554827.409547401       365.768221021
        20           10.000000       102.522367308       278.352824576   -554929.882289568   -554827.359922261       364.275226861
```

### Timings

```
 SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME
                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM
 CP2K                                 1  1.0    0.182    0.228 8802.156 8802.158
 qs_mol_dyn_low                       1  2.0    0.002    0.002 8793.361 8796.208
 qs_forces                           21  4.0    0.300    0.313 8731.917 8731.935
 qs_energies                         21  5.0    0.001    0.001 8299.259 8316.296
 ls_scf                              21  6.0    0.000    0.000 8140.304 8157.313
 velocity_verlet                     20  3.0    0.014    0.024 7689.100 7689.373
 ls_scf_main                         21  7.0    0.006    6.929 6214.282 6221.060
 dbcsr_multiply_internal           6268 10.6    5.425    5.502 5456.668 5499.570
 multiply_cannon                   6268 11.6    9.185    9.602 4797.967 4876.009
 dm_ls_curvy_optimization           234  7.8    0.002    0.002 3756.519 3757.137
 optimization_step                  234  8.8    0.001    0.046 3056.137 3059.618
 multiply_cannon_multrec         401152 12.6  845.768 1008.921 2235.185 2348.598
 mp_waitall_1                   3329028 13.7 2028.515 2273.110 2028.515 2273.110
 compute_direction_newton            78  9.8    0.101    0.977 2206.565 2206.643
 ls_scf_dm_to_ks                    255  7.9    0.003    0.004 2046.532 2047.273
 multiply_cannon_metrocomm3      401152 12.6    0.957    1.092 1280.957 1642.411
 commutator_symm                   1536 10.8    0.024    0.026 1393.033 1401.600
 dbcsr_mm_accdrv_process        7955916 13.9  862.352  933.296 1051.028 1136.256
 rebuild_ks_matrix                  277  9.6    0.002    0.004 1013.434 1013.460
 qs_ks_build_kohn_sham_matrix       277 10.6    0.037    0.051 1013.433 1013.458
 pao_post_scf                        21  7.0    0.000    0.000  979.076  979.097
 multiply_cannon_metrocomm1      401152 12.6    1.188    1.293  666.837  966.589
 pao_add_forces                      21  8.0    0.022    0.036  939.777  941.085
 qs_ks_update_qs_env                256  8.9    0.002    0.002  915.333  915.356
 update_p_exp                       234  9.8    0.013    0.032  849.491  852.997
 qs_rho_update_rho                  256  8.9    0.003    0.004  732.701  747.571
 calculate_rho_elec                 256  9.9  486.068  514.736  732.698  747.568
 sum_up_and_integrate               277 11.6    0.488    0.516  732.856  733.703
 integrate_v_rspace                 277 12.6  518.242  550.497  732.367  733.214
 matrix_qs_to_ls                    276  8.0    0.142    0.144  701.236  702.223
 pao_calc_outer_grad_lnv             21  9.0    0.100    0.114  657.554  660.960
 transform_matrix_orth              333  8.8    0.010    0.011  647.406  651.982
 purify_mcweeny_orth                255 10.7    0.009    0.493  635.637  639.508
 make_m2s                         12536 11.6    5.422    5.562  568.779  607.900
 make_images                      12536 12.6   19.448   20.772  560.314  599.309
 make_images_sizes                12536 13.6    0.028    0.032  428.487  491.794
 mp_alltoall_i44                  12536 14.6  428.459  491.763  428.459  491.763
 matrix_ls_to_qs                    276  8.9    0.154    0.157  441.958  469.591
 mp_alltoall_d11v                 10853 12.8  223.253  388.573  223.253  388.573
 dbcsr_new_transposed              3211 11.1    9.536   10.585  333.787  347.908
 rs_distribute_matrix               554 12.2   26.205   29.397  187.192  341.592
 dbcsr_redistribute                2470 12.3   31.904   35.233  319.945  334.210
```