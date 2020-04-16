from maptool.util.utils import box_center
def logo():
  s=[None]*6
  s[0]="                       _              _ "
  s[1]=" _ __ ___   __ _ _ __ | |_ ___   ___ | |"
  s[2]="| '_ ` _ \ / _` | '_ \| __/ _ \ / _ \| |"
  s[3]="| | | | | | (_| | |_) | || (_) | (_) | |"
  s[4]="|_| |_| |_|\__,_| .__/ \__\___/ \___/|_|"
  s[5]="                |_|                     "
  
  for i in range(len(s)):
      box_center(ch=s[i])
