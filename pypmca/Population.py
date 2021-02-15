# -*- coding: utf-8 -*-
"""
Population: A class that holds the state of a population as it evolves in time
as well as future contributions to the population:
    - history: a list containing the population state at each time step.
      The first element of the list corresponds to time = t0.
    - future: a list containing future additions to the population.
      The first element corresponds to the last time step of the history

A population can be used for tracking expectation values (reals) or
simulated data (integers).

@author: karlen
"""

from scipy import stats
from datetime import timedelta

from pypmca.Parameter import Parameter


class Population:
    """A class that keeps state information of a population at each time step:
        - population_name : string, short descriptor. Must be unique
        - initial_value : integer, float, or Parameter object
        - description: string, optional - for documentation
        - hidden: don't show in lists of populations when keeping list short
        - color: color to plot
        - show_sim: is it appropriate to show simulated data for this distribution?
        - report_noise: should simulated data include additional noise due to reporting?
        - if True, then report_noise_par is the lower limit of range [low,1] for which a uniform
          random number is drawn to select the number from today that are reported today
          the remaining will be in future report (along with some fraction of tomorrows report)
        - if True, then report_backlog_par is the lower limit of range [low,1] for which a uniform
          random number is drawn to select the number from the backlog that are reported today
          the remaining will be in a future report
        - if True, and the integer parameter report_days with positive value is provided, use that to determine which
          days of week that reports are produced. In BC, reports are zero on Sunday and included in the
          next reporting day. Bit encoded, with Monday = bit 0... Sunday = bit 6. BC: report_days = 63
        - if True, and the integer parameter report_days with negative value is provided, use that to identify the
          situation when data is provided on a weekly basis. The value specifies the day of week which constitutes
          a new week of data: Sunday = -6, Monday = -7, Tuesday = -1, Wednesday = -2...
        - if True and report_noise_weekly is True, then backlog is released on a weekly basis to model
          states that have a large variance in their weekly death reports
    """

    def __init__(self, population_name: str, initial_value, description: str = '',
                 hidden: bool = True, color: str = 'black', show_sim: bool = False, report_noise: bool = False,
                 report_noise_par: Parameter = None, report_backlog_par: Parameter = None,
                 report_days: Parameter = None, report_noise_weekly: bool = False):
        """
        Constructor

        Parameters
        ----------
        population_name : string, short descriptor
        initial_value : integer, float, or Parameter object

        Returns
        -------
        None.

        """
        if population_name.find(',') > -1:
            raise ValueError('Error in constructing ' + self.name +
                             ': name cannot contain a comma.')
        self.name = str(population_name)
        self.description = str(description)
        self.history = None
        self.model = None
        self.parameters = {}
        self.__initialization_by_parameter = False
        if not isinstance(initial_value, (float, int)) and not isinstance(initial_value, Parameter):
            raise TypeError('Error setting initial_value to population (' +
                            self.name + ') - must be a float, int, or Parameter object')
        self.initial_value = None
        self.set_initial_value(initial_value)

        self.future = []
        self.color = color
        self.hidden = hidden
        # identify those populations for which daily contributions are meaningful
        # this is set to False if a subtraction is performed on the population (see Subtractor.py)
        self.monotonic = True
        # identify those populations for which simulated data is appropriate to show 
        self.show_sim = show_sim
        # for tracking missed reports from yesterday (reporting noise)
        self.missed_yesterday = 0

        self.__report_noise = None
        self.__report_noise_par = None
        self.__report_backlog_par = None
        self.__report_days = None
        self.report_noise_weekly = report_noise_weekly
        self.set_report_noise(report_noise, report_noise_par, report_backlog_par, report_days)

    def set_model(self, model):
        # when a connector is added to a model, the populations are informed of the parent model
        # in order to get t0 and step size information
        self.model = model

    def set_report_noise(self, report_noise, report_noise_par, report_backlog_par, report_days):

        self.__report_noise = report_noise
        if report_noise and (report_noise_par is None or
                             not isinstance(report_noise_par, Parameter)):
            raise TypeError('Error setting report_noise_par in population (' +
                            self.name + ') - it must be a Parameter object')
        self.__report_noise_par = report_noise_par
        if report_noise and report_backlog_par is not None:
            if not isinstance(report_backlog_par, Parameter):
                raise TypeError('Error setting report_backlog_par in population (' +
                                self.name + ') - it must be a Parameter object')
            self.__report_backlog_par = report_backlog_par
            self.parameters['backlog_par'] = self.__report_backlog_par
        if report_noise and report_days is not None:
            if not isinstance(report_days, Parameter):
                raise TypeError('Error setting report_days in population (' +
                                self.name + ') - it must be a Parameter object')
            elif report_days.parameter_type != 'int':
                raise TypeError('Error setting report_days in population (' +
                                self.name + ') - it must be an integer Parameter object')
            else:
                self.__report_days = report_days
                self.parameters['report_days'] = self.__report_days

        if report_noise:
            self.parameters['noise_par'] = self.__report_noise_par
        else:
            if 'noise_par' in self.parameters:
                self.parameters.pop('noise_par')
            if 'backlog_par' in self.parameters:
                self.parameters.pop('backlog_par')
            if 'report_days' in self.parameters:
                self.parameters.pop('report_days')
        # in case parameter changed, update the model list of parameters
        # unfortunately population does not have link to model
        # self.model.update_lists()

    def get_report_noise(self):
        return {'report_noise':self.__report_noise,
                'report_noise_par':self.__report_noise_par,
                'report_backlog_par':self.__report_backlog_par,
                'report_days':self.__report_days}

    def set_initial_value(self, initial_value):
        if isinstance(initial_value, (float, int)):
            self.history = [initial_value]
            self.__initialization_by_parameter = False
            if 'init_value' in self.parameters:
                self.parameters.pop('init_value')
        elif isinstance(initial_value, Parameter):
            self.history = [initial_value.get_value()]
            self.__initialization_by_parameter = True
            self.parameters['init_value'] = initial_value
        self.initial_value = initial_value

    def __str__(self):
        return self.name

    def do_time_step(self, expectations=True):
        """Perform one step in time by incrementing population number from future.
        After doing so, remove the future element
        
        if we are reporting data (not expectations) and report_noise is true
        then add the missed reports from yesterday and miss some from today

        if report_days is not None, then save all reports on non-reporting days for
        the next day

        if self.report_noise_weekly is True, then release from backlog only
        on the first reporting day each week
        """
        next_value = 0
        if self.future is not None:

            if expectations or not self.__report_noise:
                if len(self.future) > 0:
                    next_value = self.future[0]
            else:
                incoming = 0
                if len(self.future) > 0:
                    incoming = self.future[0]

                # check if today is a reporting day
                t0 = self.model.t0
                days_after = timedelta(days=int((len(self.history)+0.5)*self.model.get_time_step()))
                today = t0 + days_after
                day_of_week = today.weekday()

                # are we reporting today? (new member added in ref model 2_2)
                # also check if reporting is done only once per week, and find day of week that
                # previous week of data is released. This is set by having report_days = -1 * day of week
                # that data is released. Note for Monday, that is specified as -7, Tuesday -1, ... Sunday -6
                try:
                    rd_not_set = self.__report_days is None
                    reporting = rd_not_set
                    weekly_reporting = False
                    if not rd_not_set:
                        weekly_reporting_day = -1 * self.__report_days.get_value()
                        if weekly_reporting_day > 0 and weekly_reporting_day < 8:
                            weekly_reporting = True
                            if weekly_reporting_day == 7:
                                weekly_reporting_day = 0
                except:
                    reporting = True
                    weekly_reporting = False

                if not reporting:
                    report_days = self.__report_days.get_value()
                    reporting = (report_days & 2**day_of_week) > 0
                    # find lowest on bit: amazing hack
                    first_reporting_day_of_week = (report_days & -report_days) == 2**day_of_week
                else:
                    first_reporting_day_of_week = day_of_week == 0
                if reporting:
                    # how many will be reported from today?
                    low_edge = self.__report_noise_par.get_value()
                    frac_report = stats.uniform.rvs(loc=low_edge, scale=1. - low_edge)
                    n_report = stats.binom.rvs(incoming, frac_report)

                    # how many will be reported from the backlog?
                    n_backlog = 0
                    if not getattr(self,'report_noise_weekly',False) or first_reporting_day_of_week:
                        try:
                            low_edge = self.__report_backlog_par.get_value()
                            frac_backlog = stats.uniform.rvs(loc=low_edge, scale=1. - low_edge)
                            n_backlog = stats.binom.rvs(self.missed_yesterday, frac_backlog)
                        except:
                            n_backlog = self.missed_yesterday

                    next_value = n_report + n_backlog

                    self.missed_yesterday -= n_backlog
                    self.missed_yesterday += incoming - n_report
                else:
                    self.missed_yesterday += incoming
                    next_value = 0

        self.history.append(self.history[-1] + next_value)
        # don't let the population go negative
        if self.history[-1] < 0:
            self.history[-1] = 0

        # if data is only reported on a weekly basis and today is first day of new week,
        # spread the previous 7 days of data equally amongst the days (note: already appended new data
        # for this new week - leave that one as is). Remainder from division of 7 spread over first days
        if not expectations and self.__report_noise and weekly_reporting:
            t0 = self.model.t0
            days_after = timedelta(days=int((len(self.history) + 0.5) * self.model.get_time_step()))
            today = t0 + days_after
            day_of_week = today.weekday()
            if weekly_reporting_day == day_of_week:
                if len(self.history) > 8:
                    week_sum = self.history[-2] - self.history[-9]
                    daily_xtra = week_sum % 7
                    daily = int((week_sum - daily_xtra) / 7)
                    cumul = self.history[-9]
                    for i in range(7):
                        iloc = i - 8
                        cumul += daily
                        if i < daily_xtra:
                            cumul += 1
                        self.history[iloc] = cumul

        # remove the future element
        if self.future is not None:
            if len(self.future) > 0:
                del self.future[0]

    def reset(self):
        """Remove history and future, reinitialize history
        """
        self.future = []
        if isinstance(self.initial_value, Parameter):
            self.history = [self.initial_value.get_value()]
        else:
            self.history = [self.initial_value]

        self.missed_yesterday = 0

    def remove_history(self):
        """Replace history with array of length 1
        """
        if len(self.history) > 0:
            current_value = self.history[-1]
            self.history = [current_value]

    def scale_history(self, scale, expectations=True):
        if len(self.history) > 0:
            for i in range(len(self.history)):
                if expectations:
                    self.history[i] *= scale
                else:
                    nu = self.history[i] * scale
                    self.history[i] = int(round(nu))

    def scale_future(self, scale, expectations=True):
        if len(self.future) > 0:
            for i in range(len(self.future)):
                if expectations:
                    self.future[i] *= scale
                else:
                    nu = self.future[i] * scale
                    self.future[i] = int(round(nu))

    def update_future_expectation(self, scale, delay):
        """Include future expections, growing the future array if needed
        """
        for i in range(len(delay.future_expectations)):
            dfe = delay.future_expectations[i] * scale
            if len(self.future) > i:
                self.future[i] += dfe
            else:
                self.future.append(dfe)

    def update_future_data(self, scale, delay):
        """Include future data, growing the future array if needed
        """
        dist = stats.multinomial.rvs(int(scale), delay.future_expectations)
        i = -1
        for n_future in dist:
            i += 1
            if len(self.future) > i:
                self.future[i] += n_future
            else:
                self.future.append(n_future)

    def update_future_fast(self, value):
        """Modify the immediate future
        """
        if len(self.future) > 0:
            self.future[0] += value
        else:
            self.future.append(value)
