def load_complement(clusters_df, Qcomplement):
  from pandas import read_pickle
  newclusters_df = read_pickle(Qcomplement)
  clusters_df = clusters_df.merge(newclusters_df, on = [('info', 'file_basename')], how = "outer", indicator = True).loc[lambda x: x['_merge'] == 'left_only'].drop('_merge', axis=1)
  return clusters_df
