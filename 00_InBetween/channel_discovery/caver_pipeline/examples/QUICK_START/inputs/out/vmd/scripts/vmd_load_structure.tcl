#set dir "C:/Downloads/caver_3.0(8)/caver_3.0/examples/QUICK_START/md_snapshots"

mol load pdb ../data/6.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

