set dir "C:/Downloads/caver_3.0(8)/caver_3.0/examples/QUICK_START/md_snapshots"

mol load pdb ${dir}/1.pdb
animate read pdb ${dir}/2.pdb
animate read pdb ${dir}/3.pdb
animate read pdb ${dir}/4.pdb
animate read pdb ${dir}/5.pdb
animate read pdb ${dir}/6.pdb
animate read pdb ${dir}/7.pdb
animate read pdb ${dir}/8.pdb
animate read pdb ${dir}/9.pdb
animate read pdb ${dir}/10.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

