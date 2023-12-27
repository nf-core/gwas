process illumina2lgen {
    input:
    set file(report), file(array)

    output:
    tuple file("${output}.ped"), file("${output}.map"), emit: ped_ch
    script:
        def samplesize = params.samplesize
	def idpat      = params.idpat
        def output = report.baseName
        """
        hostname
        topbottom.py $array $report $samplesize '$idpat' $output
        """
 }
