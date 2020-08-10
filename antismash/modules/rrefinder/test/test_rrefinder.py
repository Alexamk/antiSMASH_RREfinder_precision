
import unittest

from antismash.common.hmmer import HmmerResults

from antismash.modules.rrefinder.rrefinder import is_ripp
from antismash.modules.rrefinder.html_output import will_handle


class TestRREFinder(unittest.TestCase):
    
    def setUp(self):
        self.ripps = ['bacteriocin','cyanobactin','lanthipeptide',
                     'lassopeptide','linaridin','thiopeptide','sactipeptide',
                      'proteusin','glycocin','bottromycin','microcin']
        # Make fake HmmResults objects
        self.hmm_results = []
        
        
        
    
    def test_is_ripp(self):
        for ripp in self.ripps:
            assert is_ripp(ripp)
            assert not is_ripp(ripp[1:])
            
    def test_will_handle(self):
        assert will_handle(self.ripps)
        assert not will_handle([ripp[1:] for ripp in self.ripps])

