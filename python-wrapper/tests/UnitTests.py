#!/usr/bin/env python3

import unittest, sys
sys.path.append("../../lib/python")
import varikn


class TestInterpolation(unittest.TestCase):
    
    def setUp(self):
        self.lm1 = "../../tests/data/a.arpa"
        self.lm2 = "../../tests/data/b.arpa"
        self.data = "../../tests/data/abx20.txt"

    def test_interpolation(self):
        tg = varikn.InterTreeGram([self.lm1, self.lm2], [0.6, 0.4])
        lm = varikn.Perplexity(tg, "", "", "", "")
        lm.set_init_hist(0)
        lpsum = 0.0
        with open(self.data) as fobj:
            for line in fobj:
                for token in line.split():
                    lpsum += lm.token_logprob(token)
                lm.clear_history()
        print(lpsum)

if __name__ == '__main__':
    unittest.main()
