   # Beau Duval
# STAT 369
# 2019/10/03
# BeauVariogram.py
# This file is for HW 04 a,b.

class BeauVariogram:
    """BeauVariogram,
    This class has the capability to calculate a semivariogram, select different
    types of bins, run different models, and plot the results.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    def __init__(self,x_cords,y_cords,data,distance_type = 'euclidean'):
        """__init__, initialize the variogram objects.

        The x_cords, y_cords and data are all assumed to be associated with one
        another, and hence be the same size.

        X_mean and Y_mean are calculted in either: equalDataPoints,specifiedBins,
        or binEdgesFromWidth, and are the average distance and semivarnace
        associated with the calculated bins. These functions can be ran
        interchangably to reset the X_mean and Y_mean values. X_mean and Y_mean
        can be plotted in the plotMean function.


        Arguments:
            x_cords (list): X axis Cordinates associated with data.
            y_cords (list): Y axis Cordinates associated with data.
            data (list): Variable in which to calculate semivariance with.

        Member Feilds:
            all_points (list(list)): All x and y cordinates each in their own list.
            X_distance_data (list): 1d list of all point pair distances sorted.
            Y_semivariance_data (list): 1d list of semivarance in the order of distance.
            X_mean (list): Average point along the x axis for each bin in order.
            Y_mean (list): Average point along the y axis for each bin in order.

        """
        self.all_points = self.np.array([[x,y] for x,y in zip(x_cords,y_cords)])
        self.data = self.np.array(data).reshape((len(data),1))
        if distance_type == 'euclidean':
            self.distance_function = self.euclideanDistance
        elif distance_type == 'haversine':
            self.distance_function = self.haversineDistance



        # calculate all distances, and the individual semivarance
        self.X_distance_data,self.Y_semivariance_data = self.calcPointPairs()
        self.num_pairs = self.X_distance_data.shape[0]

        self.distance_from_all_points = self.calculated_distance_pairs()

        # we will hold universall values for the variogram, but will reset them
        # in function calls.
        self.X_mean = []
        self.Y_mean = []

        # This is the h (distance) vector that will be used to fit a model.
        self.h = self.np.linspace(min(self.X_distance_data),max(self.X_distance_data),100)

    def haversineDistance(self,point1,point2):
        #https://www.movable-type.co.uk/scripts/latlong.html
        # R is earth radius in mean km
        R = 6371.0088
        # We assume that the first index is latitude, while the second is longitude
        cord1 = self.np.radians(point1)
        cord2 = self.np.radians(point2)
        lat1 = cord1[0]
        lat2 = cord2[0]
        cord2 - cord1
        rad_diff = cord2 - cord1

        # this handles when an array and a single point are inputed.
        if len(cord1.shape) > 1 or len(cord2.shape) > 1:
            if len(cord1.shape) > len(cord2.shape):
                multi_points = cord1
                lat1 = cord2[0]
            else:
                lat1 = cord1[0]
                multi_points = cord2

        a = self.np.sin(rad_diff[:,0]/2)**2 + self.np.cos(lat1)* self.np.cos(multi_points[:,0]) * self.np.sin(rad_diff[:,1]/2)**2


        c = 2 * self.np.arctan2(self.np.sqrt(a),self.np.sqrt((1-a)))
        d = R  * c
        return d

    def euclideanDistance(self,x,y):
        """euclideanDistance, calculate the euclidean distance between points.

        This function will calculate the distance between two points, a single
        point and a list of points, and two lists of points.

        Argument:
            x (list): single point or multple points to use.
            y (list): single point or multple points to use.

        Returns:
            (np.darray): Calculated distance between input points.
        """


        return self.np.sqrt(self.np.sum((x -y)**2,axis = 1))

    def calculated_distance_pairs(self):
        distance = []
        for point in self.all_points:
            distance.extend(self.distance_function(point,self.all_points))


        return self.np.array(distance)
    def calcPointPairs(self):
        """calcPointPairs, calculate the distance from each possible pair of points.

        The distance pairs for each point will be calculated, and then the given
        data will be sorted in order of distance.

        The data and all points are expected to be the same lenght.

        #TODO: Implement different function calculation.
        The number of pairs of points should be:
        ((len(x) - 1) * len(x))/2
        """
        X_distance_data = []
        Y_semivariance_data = []
        for index,point in enumerate(self.all_points[:-1]):
            # calc distance

            X_distance_data.extend(self.distance_function(self.all_points[index + 1:],point))
            # X_distance_data.extend((self.np.sqrt(self.np.sum((self.all_points[index + 1:]
            #                                          - point)**2,axis = 1))))
            # semivariance
            # TODO: Put this into a funtion
            Y_semivariance_data.extend((((self.data[index] - self.data[index + 1:])**2)
                                        *(1/2)).flatten())


        X_distance_data = self.np.array(X_distance_data)
        Y_semivariance_data = self.np.array(Y_semivariance_data)
        # sort by distance
        distance_order = self.np.argsort(X_distance_data)

        # return sorted tuple
        return X_distance_data[distance_order],Y_semivariance_data[distance_order]

    # Equal number of data points in each bin
    # part 1

    def equalDataPoints(self,num_bins = None,num_values = None):
        """equalDataPoints, calculates bins with the same number of values in them.

        Given a value of num_bins, that many bins wll be created with equal split
        of values. Given a value of num_values, that may values will be in each bin.
        Both can not be used.

        The average x value and y value within each bin are calculated and returned,
        to be used in a semivariogram. Typically for a semvariogram X_data is
        the distance between each point. Y_data is typically individual semivariance
        values. However this function should work for any set of x and y data.

        The data is assumed to be of equal length, and ordered by X.

        Arguments:
            num_bins (int): Desired number of bins, to split values among.
            num_values (int): Desired number of values to be in each bin.
        """
        self.X_mean = []
        self.Y_mean = []

        previous_value = 0
        if num_bins:
            split = int(len(self.X_distance_data)/num_bins)
        elif num_values:
            split = num_values
        else:
            raise ParameterNotSet


        for index in range(1,len(self.X_distance_data)):
            if index % split == 0:
                self.X_mean.append(self.np.mean(self.X_distance_data[previous_value
                                                                :index]))
                self.Y_mean.append(self.np.mean(self.Y_semivariance_data[previous_value:index]))
                previous_value = index

        # This handles the last bin where it might not be possible to have the same
        # number of values given the number of desired bins.
        if previous_value != index:
            self.X_mean.append(self.np.mean(self.X_distance_data[previous_value:index]))
            self.Y_mean.append(self.np.mean(self.Y_semivariance_data[previous_value:index]))

    def specifiedBins(self,bin_edges = 40):
        """specifiedBins, calculates the average of X and Y within each bin.

        Given either a list of bin edges or the number of desired bins, this function
        will calculate the average values within each bin. This has been designed
        to createa  semivariogram.

        The data is assumed to be of equal length, and ordered by X.

        Arguments:
            bin_edges (list/int): If an int is given numpy will calculate the edges,
                else the given list of edges will be used.
        """
        if type(bin_edges) == int:
            #calculate the actual edge values
            bin_edges = self.np.histogram_bin_edges(self.X_distance_data,bins = bin_edges)


        self.X_mean = []
        self.Y_mean = []
        self.up_error = []
        self.down_error = []

        for index in range(len(bin_edges)):
            # error throws, when out of index bounds, this is the last bin!
            try:
                index_in_range = self.np.where((self.X_distance_data >= bin_edges[index])&
                                          (self.X_distance_data < bin_edges[index + 1]))
            except IndexError:
                index_in_range = self.np.where((self.X_distance_data >= bin_edges[index]))
            semi_mean = self.np.mean(self.Y_semivariance_data[index_in_range])
            temp = 1.96*self.np.std(self.Y_semivariance_data[index_in_range])/self.np.sqrt(len(self.Y_semivariance_data[index_in_range]))
            self.up_error.append(semi_mean + temp)
            self.down_error.append(semi_mean - temp)

            self.Y_mean.append(semi_mean)
            self.X_mean.append(self.np.mean(self.X_distance_data[index_in_range]))

    def visuaSelectBins(self,num_clicks):
        """visuaSelectBins,click on the raw data plot to select bins.

        The raw data is plotted, then a user is allowed to click where the bin
        edges should be. The x values from the points clicked are used to create
        the bin edges.These values are sorted so that if one clicks out of order,
        the code will not crash when binning the data.

        This data will be sent to specifiedBins to put the data into the bins.

        Argument:
            num_clicks (int): Number of desired bin edges to be selected using,
                matplotlib's ginput.
        """
        #self.plt.figure(figsize=(10,7))
        self.plotRawScatter()
        clicked_points = self.plt.ginput(num_clicks)
        x,y = zip(*clicked_points)
        self.specifiedBins(bin_edges = sorted(list(x)))

    def binEdgesFromWidth(self,bin_width):
        """binEdgesFromWidth, calculates binning edges from a given bin width.

        Arguments:
            X_data (list): The data in which the bins will be created on.
            bin_width (int): Desired uniform width between bin edges in X_data units.
        """

        max_distance = max(self.X_distance_data)
        min_distance = min(self.X_distance_data)

        num_bins = int((max_distance - min_distance)/bin_width)
        bin_edges = self.np.linspace(min_distance,max_distance,num_bins + 1)

        # now that we know the bin edges, select the data that should go in them
        self.specifiedBins(bin_edges)

    def plotRawScatter(self,opacity = 1):
        """plotRawScatter, will plot the raw distance and semivarance data.

        These data points are initially calculated once the variogram object is
        initiated. See function calcPointPairs.

        Argument:
            opacity (float): Alpha value to set that points on the plot to.
        """
        self.plt.scatter(self.X_distance_data,self.Y_semivariance_data,alpha = opacity)


    def plotMean(self):
        """plotMean, will plot the mean x and y values within the calculated bins.
        """
        self.plt.scatter(self.X_mean,self.Y_mean)


    def setModelParams(self,sill,nugget,a):
        """setModelParams, define the paramaters for all semivariogram Models.
        This was done to allow for multiple models to be easily ran without
        setting the same model paramaters over and over agian.

        These paramaters are defined in Introduction to Applied Geostatistic
        page 143 chapter 7.

        !Warning! the term range was not used here in the code, because it is a
        key word in python, and will disrupt other code.

        Arguments:
            sill (float): The plateau the variogram reaches at the range.
            nugget (float): The vertical jump from the origin to the first
                values on the variogram.
            a (float): The distanc at which the variogram plateaus, aslso known
                as the range.
        """
        self.sill = sill
        self.nugget = nugget
        self.a = a
    def setModel(self,model_function):
        self.model_function = model_function

    def modelTransform(self,data):
        return self.model_function(data)

    def exponential(self,transform = False):
        """exponential, semivariogram model as defined in STATS 369.
        See de Marsily for refrence.

        This function calculates the model, and plots it with the according label.
        To be used with a new figure.

        Nugget, Sill, and a (range) are all defined in setModelParams.

        h is defined as a vector along the x axis in init, in which to calculate
        the model on.
        """
        # This is to be used to transform a specific vector of values.
        if isinstance(transform,type(self.np.array([]))):
            return self.nugget +(self.sill-self.nugget)*(1 - self.np.exp((-3*transform)/self.a))
        # This is used for plotting.
        else:
            y_h = self.nugget +(self.sill-self.nugget)*(1 - self.np.exp((-3*self.h)/self.a))
            self.plt.plot(self.h,y_h,label = "Exponential Model")

    def gaussian(self):
        """gaussian, semivariogram model as defined in STATS 369.
        See de Marsily for refrence.

        This function calculates the model, and plots it with the according label.
        To be used with a new figure.

        Nugget, Sill, and a (range) are all defined in setModelParams.

        h is defined as a vector along the x axis in init, in which to calculate
        the model on.
        """
        y_h = self.nugget +(self.sill-self.nugget)*(1 - self.np.exp(-((self.np.sqrt(3)*self.h)/self.a)**2))
        self.plt.plot(self.h,y_h,label = 'Gaussian Model')

    def spherical(self,transform = False):
        """spherical, semivariogram model as defined in STATS 369.
        See de Marsily for refrence.

        This function calculates the model, and plots it with the according label.
        To be used with a new figure.

        h is defined as a vector along the x axis in init, in which to calculate
        the model on.

        Nugget, Sill, and a (range) are all defined in setModelParams.
        """
        if isinstance(transform,type(self.np.array([]))):
            temp = transform[transform <= self.a]
            y_h = self.nugget +(self.sill-self.nugget)*(1.5*(temp/self.a)-0.5*(temp/self.a)**3)
            y_h = self.np.concatenate((y_h,self.np.full(transform[transform > self.a].shape,
                                                        self.nugget +(self.sill-self.nugget))))
            return y_h
        else:
            temp = self.h[self.h <= self.a]
            y_h = self.nugget +(self.sill-self.nugget)*(1.5*(temp/self.a)-0.5*(temp/self.a)**3)
            y_h = self.np.concatenate((y_h,self.np.full(self.h[self.h > self.a].shape,
                                                        self.nugget +(self.sill-self.nugget))))
            self.plt.plot(self.h,y_h,label = 'Spherical Model')

    def linear(self):
        """linaer, semivariogram model as defined in STATS 369.
        See de Marsily for refrence.

        This function calculates the model, and plots it with the according label.
        To be used with a new figure.

        h is defined as a vector along the x axis in init, in which to calculate
        the model on.

        Nugget, Sill, and a (range) are all defined in setModelParams.
        """
        y_h = self.nugget +(self.sill-self.nugget)*(self.h/self.a)
        self.plt.plot(self.h,y_h,label = 'Linaer Model')

    def cubic(self):
        """cubic, semivariogram model as defined in STATS 369.
        See de Marsily for refrence.

        This function calculates the model, and plots it with the according label.
        To be used with a new figure.

        h is defined as a vector along the x axis in init, in which to calculate
        the model on.

        Nugget, Sill, and a (range) are all defined in setModelParams.
        """

        y_h = self.nugget +(self.sill-self.nugget) * (7*(self.h/self.a)**2 - 8.75*(self.h/self.a)**3
                                       + 3.5*(self.h/self.a)**5 - 0.75*(self.h/self.a)**7)
        self.plt.plot(self.h,y_h,label = 'Cubic Model')

    def showFig(self,legend = False,save_fig_name = None,xlabel_text = None):
        """showFig, will display the figure created from newFig.

        figure axi labels are added on in here.
        showFig also has the capability to allow for a legend, and saving the
        figure

        Arguments:
            legend (bool): If one wants a legend, should be used when fitting models.
            save_fig_name (str): Name of file that will be saved to Plots folder.
        """
        self.plt.xlabel("{}".format(xlabel_text)) if xlabel_text else self.plt.xlabel("Distance")
        self.plt.ylabel("Semivariance")
        self.plt.legend() if legend else None
        self.plt.savefig("Plots/{}.png".format(save_fig_name)) if save_fig_name else None
        self.plt.show(self.figure)

    def newFig(self,title_text = None,dpi = None):
        """newFig, creates a new figure object for internal use.

        This should be used each time before plotting. However data will plot
        and be shown even if this is not used, it will just be formated poorly.

        Argument:
            title_text (str): Desired Title for new figure.
        """
        if dpi:
            self.figure = self.plt.figure(figsize = (10,7),dpi = dpi)
        else:
            self.figure = self.plt.figure(figsize = (10,7))
        self.plt.title(title_text) if title_text else None

class ParameterNotSet(Exception):
    """ParameterNotSet, is for internal use in equalDataPoints.

    raise an error and show specific message.
    """
    def __init__(self):
        Exception.__init__(self,'Either num_bins or num_values must be set.')
