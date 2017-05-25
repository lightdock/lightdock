
class GSOClusterTask(object):
    """A GSO execution in a given cluster"""
    def __init__(self, id_cluster, gso, steps, dest_folder):
        self.id = id_cluster
        self.gso = gso
        self.steps = steps
        self.saving_path = dest_folder

    def run(self):
        self.gso.run(self.steps, cluster_id=self.id, verbose=True,
                     saving_path=self.saving_path, save_intermediary=True)
