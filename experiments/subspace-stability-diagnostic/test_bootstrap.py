"""Standard-library packaging tests; do not start a numerical worker."""
from pathlib import Path
import json,os,sys,tempfile,unittest
import bootstrap_replay as b
class BootstrapTests(unittest.TestCase):
 def test_portable_fresh_protocol_and_history_preservation(self):
  before=b.digest(b.ROOT/'replay-template.json')
  with tempfile.TemporaryDirectory(dir=b.ROOT,prefix='.bootstrap-test-') as tmp:
   out=Path(tmp)/'fresh';receipt=b.prepare(out,sys.executable,sys.executable);p=json.loads((out/'protocol.json').read_text())
   self.assertEqual(receipt['status'],'PREPARED; NOT EXECUTED');self.assertEqual(len(p['jobs']),23);self.assertFalse((out/'execution.json').exists());self.assertEqual(p['executables']['pinned'],os.path.abspath(sys.executable));self.assertEqual(p['prior_execution_seconds'],0)
   self.assertEqual((out/'protocol.sha256').read_text().strip(),b.digest(out/'protocol.json'))
   for name,expected in p['code_sha256'].items():self.assertEqual(b.digest(out/name),expected)
   with self.assertRaises(FileExistsError):b.prepare(out,sys.executable,sys.executable)
  self.assertEqual(b.digest(b.ROOT/'replay-template.json'),before)
 def test_missing_interpreter_does_not_create_output(self):
  with tempfile.TemporaryDirectory(dir=b.ROOT,prefix='.bootstrap-test-') as tmp:
   out=Path(tmp)/'fresh'
   with self.assertRaises(FileNotFoundError):b.prepare(out,Path(tmp)/'missing',sys.executable)
   self.assertFalse(out.exists())
 def test_virtual_environment_symlink_path_is_preserved(self):
  with tempfile.TemporaryDirectory(dir=b.ROOT,prefix='.bootstrap-test-') as tmp:
   link=Path(tmp)/'venv/bin/python';link.parent.mkdir(parents=True);link.symlink_to(Path(sys.executable).resolve());out=Path(tmp)/'fresh'
   b.prepare(out,str(link),sys.executable);p=json.loads((out/'protocol.json').read_text())
   self.assertEqual(p['executables']['pinned'],str(link.absolute()))
   self.assertNotEqual(p['executables']['pinned'],str(link.resolve()))
if __name__=='__main__':unittest.main()
