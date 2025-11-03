import pandas as pd

class ResultAggregator:
    def __init__(self):
        self.results = []
    
    def add_result(self, result):
        self.results.append(result)
    
    def to_dataframe(self):
        return pd.DataFrame(self.results)
    
    def save_csv(self, filepath):
        df = self.to_dataframe()
        df.to_csv(filepath, index=False)
        print(f"Results saved to {filepath}")
