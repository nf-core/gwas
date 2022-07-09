def getres(x) {
  def command1 = "$x"
  def command2 = "head -n 1"
  def proc1 = command1.execute()
  def proc2 = command2.execute()
  def proc = proc1 | proc2
  proc.waitFor()
  def res ="${proc.in.text}"
  return res.trim()
}


def checkColumnHeader(fname, columns) {
  if (workflow.profile == "awsbatch") return;
  if (fname.toString().contains("s3://")) return;
  if (nullfile.contains(fname)) return;
  new File(fname).withReader { line = it.readLine().tokenize() }
  problem = false;
  columns.each { col ->
    if (! line.contains(col) ) {
      println "The file <$fname> does not contain the column <$col>";
      problem=true;
    }
    if (problem)
      System.exit(2)
  }
}



def checkSampleSheet(fname)  {
  if (workflow.profile == "awsbatch") return;
  if (fname.contains("s3://") )return;
  if (fname.contains("az://") ) return;
  if (nullfile.contains(fname) || fname.contains(".xls")) return;
  new File(fname).withReader { line = it.readLine()}
  problem  = false
  prob_str = ""
  if (! line.contains(",")) {
    problem = true;
    prob_str = "If given as a CSV file, it must be comma-separated\n";
  }
  headers = line.tokenize(",")
  headers.each { println it}
  if (!(headers.contains("Institute Sample Label") ||
      (headers.contains("Sample Plate") && headers.contains("Well")))) {
    problem= true
    prob_str = prob_str + "Column headers must include 'Institute Sample Label'  or both 'Sample Plate' and 'Well'"
  }
  if (problem)  {
    println "There's a problem with the sample sheet <$fname>."
    println prob_str;
    //FIXME
    System.exit(1)
  }
}




// This method first checks that the data file has the stated column
// If so, it creates a channel for it
// NB: if the file is in S3 we cannot do the test since Groovy does not
// allow us to access the file directly
def getSubChannel = { parm, parm_name, col_name ->
  if (parm.toString().contains("s3://")) {
    println "The file <$parm> is in S3 so we cannot do a pre-check";
    return Channel.fromPath(parm);
  }
  if (parm.toString().contains("az://")) {
    println "The file <$parm> is in Azure so we cannot do a pre-check";
    return Channel.fromPath(parm);
  }
  if ((parm==0) || (parm=="0") || (parm==false) || (parm=="false")) {
    filename = "emptyZ0${parm_name}.txt";
    new File(filename).createNewFile()
    new_ch = Channel.fromPath(filename);

  } else {
    if (! file(parm).exists()) {
     error("\n\nThe file <$parm> given for <params.${parm_name}> does not exist")
    } else {
      def line
      new File(parm).withReader { line = it.readLine() }
      fields = line.split()
      if (! fields.contains(col_name))
	  error("\n\nThe file <$parm> given for <params.${parm_name}> does not have a column <${col_name}>\n")
    }
    new_ch = Channel.fromPath(parm);
  }
  return new_ch;
}
