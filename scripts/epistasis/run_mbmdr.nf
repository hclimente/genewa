ped = file("$params.ped")
map = file("$params.map")

N = 10

process addMapHeader {

	input:
		file map

	output:
		file "${map.baseName}.fixed.map" into fixedMap

	"""
	echo -e 'chr\tname\tcM\tpos' >${map.baseName}.fixed.map
	cat $map >>${map.baseName}.fixed.map
	"""

}

process ped2mbmdr {

	input:
		file ped
		file fixedMap

	output:
		file "${ped.baseName}.mdr" into mbmdrIn
		file "${ped.baseName}.labels" into mbmdrLabels

	"""
	mbmdr.out --plink2mbmdr --binary -ped $ped -map $fixedMap -o ${ped.baseName}.mdr -tr ${ped.baseName}.labels
	"""

}

mbmdrIn.into { mbmdrIn_1; mbmdrIn_2; mbmdrIn_3; mbmdrIn_4 }

process run_mbmdr_1 {

	input:
		each i from 1..N
		file mbmdrIn_1

	output:
		file "top*.txt" into topFiles

	"""
	mbmdr.out --binary --gammastep1 -i $i -N $N $mbmdrIn_1
	"""

}

process run_mbmdr_2 {

	input:
		file "top*.txt" from topFiles.collect()
		file mbmdrIn_2

	output:
		file "topFile.txt" into topFile

	"""
	mbmdr.out --binary --gammastep2 -N $N $mbmdrIn_2
	"""

}

topFile.into { topFile_3; topFile_4 }

process run_mbmdr_3 {

	input:
		each i from 1..N
		file topFile_3
		file mbmdrIn_3

	output:
		file "perm_*" into permutFiles

	"""
	mbmdr.out --binary --gammastep3 -q $i -o perm_$i -t $topFile_3 $mbmdrIn_3
	"""

}

process run_mbmdr_4 {

	input:
		file "perm_*.txt" from permutFiles.collect()
		file topFile_4
		file mbmdrIn_4

	output:
		file "*_output.txt" into output

	"""
	mbmdr.out --binary --gammastep4 -c perm_ -q ${permutFiles.toList().size()} -t $topFile_4 $mbmdrIn_4
	"""

}
