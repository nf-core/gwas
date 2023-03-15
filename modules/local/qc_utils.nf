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
