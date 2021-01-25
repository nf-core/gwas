params.vcfFile = null

rawFileChannel = Channel.fromPath(params.vcfFile)

process check_vcf_build {

  echo true

  input:
  file vcfFile from rawFileChannel
  output:
  file("*.BuildChecked") into inputFileOutChan

  script:
    """
    /app/0_check_vcf_build.slurm.sh $vcfFile .
    """
}
