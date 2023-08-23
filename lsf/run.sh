#!/bin/bash

bsub -W 180 -n 1 Rscript results_abiotic.r
for species in 'Amerikaanse eik' 'Asperge_Liggende asperge' 'Bevertjes' 'Bochtige smele' 'Bosaardbei' 'Bosanemoon' 'Brede stekelvaren' 'Brede wespenorchis_Duinwespenorchis' 'Dalkruid' 'Duizendblad' 'Eenarig wollegras' 'Engels raaigras' 'Es' 'Gestreepte witbol' 'Gevlekte orchis_Bosorchis' 'Gewone dophei' 'Gewone vogelkers' 'Gewoon kweldergras' 'Gewoon reukgras' 'Grote brandnetel' 'Hazelaar' 'Hulst' 'Jakobskruiskruid_Duinkruiskruid' 'Kleine ratelaar' 'Kleine veenbes' 'Kleine zonnedauw' 'Knolboterbloem' 'Kraaihei' 'Kruipwilg' 'Moeraswespenorchis' 'Muurpeper' 'Parnassia' 'Peen' 'Pitrus' 'Ratelpopulier' 'Ruwe berk' 'Ruwe smele' 'Slanke sleutelbloem' 'Smalle stekelvaren' 'Spaanse aak' 'Tengere rus' 'Tweerijige zegge' 'Veenpluis' 'Watermunt' 'Wilde gagel' 'Witte klaverzuring' 'Zachte berk' 'Zandstruisgras' 'Zoete kers' 'Zomereik'
do
  bsub -W 1440 -n 1 Rscript results_certain.r "$species"
  bsub -W 1440 -n 1 Rscript results_uncertain.r "$species"
  sleep 1
done

for species in 'Amerikaanse eik' 'Asperge_Liggende asperge' 'Bevertjes' 'Bochtige smele' 'Bosaardbei' 'Bosanemoon' 'Brede stekelvaren' 'Brede wespenorchis_Duinwespenorchis' 'Dalkruid' 'Duizendblad' 'Eenarig wollegras' 'Engels raaigras' 'Es' 'Gestreepte witbol' 'Gevlekte orchis_Bosorchis' 'Gewone dophei' 'Gewone vogelkers' 'Gewoon kweldergras' 'Gewoon reukgras' 'Grote brandnetel' 'Hazelaar' 'Hulst' 'Jakobskruiskruid_Duinkruiskruid' 'Kleine ratelaar' 'Kleine veenbes' 'Kleine zonnedauw' 'Knolboterbloem' 'Kraaihei' 'Kruipwilg' 'Moeraswespenorchis' 'Muurpeper' 'Parnassia' 'Peen' 'Pitrus' 'Ratelpopulier' 'Ruwe berk' 'Ruwe smele' 'Slanke sleutelbloem' 'Smalle stekelvaren' 'Spaanse aak' 'Tengere rus' 'Tweerijige zegge' 'Veenpluis' 'Watermunt' 'Wilde gagel' 'Witte klaverzuring' 'Zachte berk' 'Zandstruisgras' 'Zoete kers' 'Zomereik'
do
  for provincie in 'Drenthe' 'Fryslân' 'Overijssel' 'Flevoland' 'Gelderland' 'Utrecht' 'Groningen' 'Limburg' 'Noord-Brabant' 'Zuid-Holland' 'Zeeland' 'Noord-Holland'
  do
    bsub -W 1440 -n 1 Rscript predict_certain.r "$species" "$provincie"
    bsub -W 1440 -n 1 Rscript predict_uncertain.r "$species" "$provincie"
    sleep 1
  done
done

for species in 'Amerikaanse eik' 'Asperge_Liggende asperge' 'Bevertjes' 'Bochtige smele' 'Bosaardbei' 'Bosanemoon' 'Brede stekelvaren' 'Brede wespenorchis_Duinwespenorchis' 'Dalkruid' 'Duizendblad' 'Eenarig wollegras' 'Engels raaigras' 'Es' 'Gestreepte witbol' 'Gevlekte orchis_Bosorchis' 'Gewone dophei' 'Gewone vogelkers' 'Gewoon kweldergras' 'Gewoon reukgras' 'Grote brandnetel' 'Hazelaar' 'Hulst' 'Jakobskruiskruid_Duinkruiskruid' 'Kleine ratelaar' 'Kleine veenbes' 'Kleine zonnedauw' 'Knolboterbloem' 'Kraaihei' 'Kruipwilg' 'Moeraswespenorchis' 'Muurpeper' 'Parnassia' 'Peen' 'Pitrus' 'Ratelpopulier' 'Ruwe berk' 'Ruwe smele' 'Slanke sleutelbloem' 'Smalle stekelvaren' 'Spaanse aak' 'Tengere rus' 'Tweerijige zegge' 'Veenpluis' 'Watermunt' 'Wilde gagel' 'Witte klaverzuring' 'Zachte berk' 'Zandstruisgras' 'Zoete kers' 'Zomereik'
do
  for geo in 'afz' 'duo' 'duw' 'hll' 'hzn' 'hzo' 'hzv' 'hzz' 'lvh' 'lvn' 'riv' 'zkm' 'zkn' 'zkz'
  do
    bsub -W 1440 -n 1 Rscript predict_certain2.r "$species" "$geo"
    bsub -W 1440 -n 1 Rscript predict_uncertain2.r "$species" "$geo"
    sleep 1
  done
done

# After the jobs are finished...
# cat /mnt/scratch_dir/viljanem/results_certain/*.csv /mnt/scratch_dir/viljanem/results_uncertain/*.csv > results.csv
# cat /mnt/scratch_dir/viljanem/predictions_certain/*.csv /mnt/scratch_dir/viljanem/predictions_uncertain/*.csv > predictions_province.csv 
# cat /mnt/scratch_dir/viljanem/prevalences_certain/*.csv /mnt/scratch_dir/viljanem/prevalences_uncertain/*.csv > prevalences_province.csv 
# cat /mnt/scratch_dir/viljanem/predictions_certain2/*.csv /mnt/scratch_dir/viljanem/predictions_uncertain2/*.csv > predictions_fgr.csv 
# cat /mnt/scratch_dir/viljanem/prevalences_certain2/*.csv /mnt/scratch_dir/viljanem/prevalences_uncertain2/*.csv > prevalences_fgr.csv 
# cp *.csv /data/BioGrid/viljanem/sprr/
