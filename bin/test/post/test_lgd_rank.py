"""Test for lgd_rank post script"""

import os
import filecmp
import shutil
from ..regression import RegressionTest


class TestGenerateRanking(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_lgd_rank/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/4IZ7/'

    def teardown(self):
        self.clean_test_path()

    def test_rank_with_clusters(self):
        num_swarms = 4
        num_steps = 10

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        for i in range(num_swarms):
            swarm_dir = 'swarm_%d' % i
            os.mkdir(swarm_dir)
            shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'gso_%d.out' % num_steps),
                            os.path.join(self.test_path, swarm_dir, 'gso_%d.out' % num_steps))
            shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'cluster.repr'),
                            os.path.join(self.test_path, swarm_dir, 'cluster.repr'))

        
        command = "lgd_rank.py %d %d > test.out" % (num_swarms, num_steps)
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'rank_by_scoring.list'), 
                           os.path.join(self.test_path, 'rank_by_scoring.list'))


class TestGenerateRankingWithoutClusters(RegressionTest):

    def setup(self):
        self.path = os.path.dirname(os.path.realpath(__file__))
        self.test_path = self.path + '/scratch_lgd_rank_no_clust/'
        self.ini_test_path()
        self.golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + \
                                '/golden_data/4IZ7/'

    def teardown(self):
        self.clean_test_path()

    def test_rank_with_clusters(self):
        num_swarms = 4
        num_steps = 10

        # Prepare folder structure for this test
        os.chdir(self.test_path)
        for i in range(num_swarms):
            swarm_dir = 'swarm_%d' % i
            os.mkdir(swarm_dir)
            shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'gso_%d.out' % num_steps),
                            os.path.join(self.test_path, swarm_dir, 'gso_%d.out' % num_steps))
            shutil.copyfile(os.path.join(self.golden_data_path, swarm_dir, 'cluster.repr'),
                            os.path.join(self.test_path, swarm_dir, 'cluster.repr'))

        
        command = "lgd_rank.py %d %d --ignore_clusters > test.out" % (num_swarms, num_steps)
        os.system(command)

        assert filecmp.cmp(os.path.join(self.golden_data_path, 'rank_by_scoring_noclust.list'), 
                           os.path.join(self.test_path, 'rank_by_scoring.list'))
