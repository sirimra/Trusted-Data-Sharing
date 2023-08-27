# Importing the geopy.distance module for distance calculations
import geopy.distance

class PointClass():
    # Class attributes to store point ID and distance using geopy.distance
    point_id : int
    distance : geopy.distance

    # Constructor to initialize point ID and distance
    def __init__(self, point_id, distance):
        self.point_id=point_id
        self.distance= distance


    # Method to convert a string expression to point attributes
    def convert_string(self, string_exp):
        res = eval (string_exp)
        self.point_id = res[0]
        self.distance =geopy.distance.distance(res[1])


    # Method to convert the PointClass instance to a string
    def to_str(self):
        return str(self.point_id)+','+str(self.distance.km)


