import requests
import json
import polyline

class RouteBuilder:
    def __init__(self):
        self.routes = []

    def makeRequest(self, coords1, coords2, returnRoute=False):
        #makes single request and stores route in a list; does not clear route list
        lat_1, lon_1 = coords1
        lat_2, lon_2 = coords2
        # call the OSMR API
        r = requests.get(f"http://router.project-osrm.org/route/v1/driving/{lon_1},{lat_1};{lon_2},{lat_2}?overview=simplified")
        # then you load the response using the json libray
        # by default you get only one alternative so you access 0-th element of the `routes`
        #routes = json.loads(r.content)
        r.json()
        #self.routes.append(r.json()['routes'][0])
        return r.json()['routes'][0]

    def makeRequest1(self, coords1, coords2, returnRoute=False):
        #makes single request and stores route in a list; does not clear route list
        lat_1, lon_1 = coords1
        lat_2, lon_2 = coords2
        # call the OSMR API
        r = requests.get(f"http://router.project-osrm.org/route/v1/driving/{lon_1},{lat_1};{lon_2},{lat_2}?overview=simplified")
        # then you load the response using the json libray
        # by default you get only one alternative so you access 0-th element of the `routes`
        #routes = json.loads(r.content)
        return r.json()['routes'][0]["distance"]

    def makeRequests(self, coord_lst):
        #multiple request; cleare route list before adding routes
        self.routes = []
        for x in coord_lst:
            self.makeRequest(x)

    def getDists(self):
        return [route["distance"] for route in self.routes]

    def getTimes(self):
        return [route["duration"] for route in self.routes]

    def getPolylines(self):
        return [polyline.decode(route['geometry']) for route in self.routes]
