# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:17:54 2017

Author: Paul TO, Ka Po
Contact: kpto@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Division of Biomedical Engineering

Supervisor: Prof. Henry LAM, H. N.
Contact: kehlam@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Department of Chemical and Biomolecular Engineering

Desciption of this module:
"""


# ====BEGIN OF MODULE IMPORT====
import logging
from ast import literal_eval
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
RANGEABLE_DATA_TYPE = (int, float)
INF = float('inf')
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class ParameterBase:
    def __init__(self, default, type_, range_=(-INF, INF), immutable=False):
        self._value_type = type_
        self._immutable = False
        if type_ not in RANGEABLE_DATA_TYPE:
            self._value_range = None
        else:
            self._value_range = range_
        self._value = None
        self.value = default
        self._immutable = immutable
        return

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, new_value):
        if self._immutable:
            err_msg = '\nThis parameter is set to be immutable.' +\
                      '\nIf you insist to change it, make a new instance of {}.'.format(type(self)) +\
                      '\nMake sure that you know what you are doing!!'
            logging.error(err_msg)
            raise TypeError(err_msg)
        raw = new_value
        if self._value_type is not str and type(new_value) is str:
            try:
                new_value = literal_eval(new_value)
            except Exception:
                err_msg = '\nInput: {}'.format(raw) +\
                          '\nUnable to parse the input, input must follow Python syntax.' +\
                          '\nFollowing objects are supported by literal_eval:' +\
                          '\nstrings, bytes, numbers, tuples, lists, dicts, sets, booleans, and None.'
                logging.error(err_msg)
                raise
        self._valid_type(new_value)
        self._valid_value(new_value)
        self._other_validation(new_value)
        self._value = new_value
        return

    def _valid_type(self, value):
        if type(value) != self._value_type:
            err_msg = '\nInvalid value type.' +\
                      '\nInput: {}'.format(type(value)) +\
                      '\nExpect: {}'.format(self._value_type)
            logging.error(err_msg)
            raise TypeError(err_msg)
        return

    def _valid_value(self, value):
        if self._value_range is not None and not with_in_range(self._value_range, value):
            err_msg = '\nInvalid value.' +\
                      '\nInput: {}'.format(value) +\
                      '\nExpect value between {} - {}'\
                      .format(self._value_range[0], self._value_range[1])
            logging.error(err_msg)
            raise ValueError(err_msg)
        return

    # virtual function, overridden by sub-class when necessary
    def _other_validation(self, value):
        return

    def is_immutable(self):
        return self._immutable

    def get_value_type(self):
        return self._value_type

    def get_value_range(self):
        if self._value_range is None:
            err_msg = '\nValue range is not applicable for this parameter with value type of {}.'\
                      .format(self._value_type)
            logging.error(err_msg)
            raise TypeError(err_msg)
        return self._value_range

    def __str__(self):
        return str(self._value)


# virtual class, instance should never be created.
class ParameterElements(ParameterBase):
    _type = None

    def __init__(self, default, element_type, element_range=(-INF, INF), immutable=False):
        self._element_type = element_type
        if element_type not in RANGEABLE_DATA_TYPE:
            self._element_range = None
        else:
            self._element_range = element_range
        super().__init__(default, type_=type(self)._type, range_=None, immutable=immutable)
        return

    def _other_validation(self, value):
        self._element_validation(value)
        return True

    def _element_validation(self, value):
        for e in value:
            # element type checking
            if type(e) != self._element_type:
                err_msg = '\nInput: {}'.format(value) +\
                          '\nCatched invalid element: {} {}'.format(e, type(e)) +\
                          '\nExpected type of elements is {}.'.format(self._element_type)
                logging.error(err_msg)
                raise ValueError(err_msg)
            # element range checking
            if self._element_range is not None and not with_in_range(self._element_range, e):
                err_msg = '\nInput: {}'.format(value) +\
                          '\nCatched invalid element: {}'.format(e) +\
                          '\nElements should be values between {} - {}.'\
                          .format(self._element_range[0], self._element_range[1])
                logging.error(err_msg)
                raise ValueError(err_msg)
        return

    def get_value_range(self):
        err_msg = '\nValue range is not applicable for {}.'.format(type(self)) +\
                  '\nYou may be looking for element_range.'
        logging.error(err_msg)
        raise TypeError(err_msg)

    def get_element_type(self):
        return self._element_type

    def get_element_range(self):
        if self._element_range is None:
            err_msg = '\nElement range is not applicable for this parameter with element type of {}.'\
                      .format(self._element_type)
            logging.error(err_msg)
            raise TypeError(err_msg)
        return self._element_range


class ParameterTuple(ParameterElements):
    _type = tuple

    def __init__(self, default, element_type, element_range=(-INF, INF), size=INF, immutable=False):
        self._size = size
        super().__init__(default, element_type=element_type, element_range=element_range, immutable=immutable)
        return

    def _other_validation(self, value):
        self._size_validation(value)
        self._element_validation(value)
        return True

    def _size_validation(self, value):
        if self._size is not INF and len(value) != self._size:
            err_msg = '\nInput: {}'.format(value) +\
                      '\nInput has the length of {}.'.format(len(value)) +\
                      '\nThis parameter must be a tuple with size of {}.'.format(self._size)
            logging.error(err_msg)
            raise ValueError(err_msg)
        return

    def get_size(self):
        return self._size


class ParameterRange(ParameterTuple):
    _type = tuple

    def __init__(self, default, element_type, element_range=(-INF, INF), immutable=False):
        super().__init__(default, element_type=element_type, element_range=element_range,
                         size=2, immutable=immutable)
        return

    def _other_validation(self, value):
        self._size_validation(value)
        self._element_validation(value)
        self._range_validation(value)
        return

    # check order of lower and upper limits
    def _range_validation(self, value):
        if value[1] < value[0]:
            err_msg = '\nInput: {}'.format(value) +\
                      '\nThe upper limit is smaller than the lower limit.' +\
                      '\nFirst element is the lower limit, second element is the upper limit.'
            logging.error(err_msg)
            raise ValueError(err_msg)


class ParameterList(ParameterElements):
    _type = list
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
def with_in_range(range_, value):
    return range_[0] <= value <= range_[1]
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
