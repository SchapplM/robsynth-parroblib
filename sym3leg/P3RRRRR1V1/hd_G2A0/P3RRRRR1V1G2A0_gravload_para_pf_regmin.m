% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR1V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x18]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR1V1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:35:00
% EndTime: 2020-08-07 03:35:01
% DurationCPUTime: 1.05s
% Computational Cost: add. (873->162), mult. (1434->279), div. (180->8), fcn. (1542->36), ass. (0->156)
t2011 = legFrame(2,2);
t1995 = sin(t2011);
t1998 = cos(t2011);
t1983 = t1998 * g(1) - t1995 * g(2);
t2018 = sin(qJ(1,2));
t2027 = cos(qJ(1,2));
t1959 = -g(3) * t2018 + t1983 * t2027;
t2008 = qJ(2,2) + qJ(3,2);
t1992 = cos(t2008);
t2026 = cos(qJ(2,2));
t1974 = 0.1e1 / (t2026 * pkin(2) + pkin(3) * t1992 + pkin(1));
t2089 = t1974 * t2027;
t2060 = t1998 * t2089;
t2047 = t1959 * t2060;
t2061 = t1995 * t2089;
t2048 = t1959 * t2061;
t2090 = t1974 * t2018;
t2070 = t1959 * t2090;
t2012 = legFrame(1,2);
t1996 = sin(t2012);
t1999 = cos(t2012);
t1984 = t1999 * g(1) - t1996 * g(2);
t2021 = sin(qJ(1,1));
t2030 = cos(qJ(1,1));
t1960 = -g(3) * t2021 + t1984 * t2030;
t2009 = qJ(2,1) + qJ(3,1);
t1993 = cos(t2009);
t2029 = cos(qJ(2,1));
t1975 = 0.1e1 / (t2029 * pkin(2) + pkin(3) * t1993 + pkin(1));
t2087 = t1975 * t2030;
t2058 = t1999 * t2087;
t2043 = t1960 * t2058;
t2059 = t1996 * t2087;
t2044 = t1960 * t2059;
t2088 = t1975 * t2021;
t2068 = t1960 * t2088;
t2010 = legFrame(3,2);
t1994 = sin(t2010);
t1997 = cos(t2010);
t1982 = t1997 * g(1) - t1994 * g(2);
t2015 = sin(qJ(1,3));
t2024 = cos(qJ(1,3));
t1961 = -g(3) * t2015 + t1982 * t2024;
t2007 = qJ(2,3) + qJ(3,3);
t1991 = cos(t2007);
t2023 = cos(qJ(2,3));
t1973 = 0.1e1 / (t2023 * pkin(2) + pkin(3) * t1991 + pkin(1));
t2092 = t1973 * t2015;
t2107 = t1961 * t2092;
t1979 = t1994 * g(1) + t1997 * g(2);
t1988 = sin(t2007);
t2053 = g(3) * t2024 + t1982 * t2015;
t1943 = t1979 * t1991 + t1988 * t2053;
t2013 = sin(qJ(3,3));
t2004 = 0.1e1 / t2013;
t2104 = t1943 * t2004;
t1944 = -t1979 * t1988 + t1991 * t2053;
t2103 = t1944 * t2004;
t1980 = t1995 * g(1) + t1998 * g(2);
t1989 = sin(t2008);
t2052 = g(3) * t2027 + t1983 * t2018;
t1945 = t1980 * t1992 + t1989 * t2052;
t2016 = sin(qJ(3,2));
t2005 = 0.1e1 / t2016;
t2102 = t1945 * t2005;
t1946 = -t1980 * t1989 + t1992 * t2052;
t2101 = t1946 * t2005;
t1981 = t1996 * g(1) + t1999 * g(2);
t1990 = sin(t2009);
t2051 = g(3) * t2030 + t1984 * t2021;
t1947 = t1981 * t1993 + t1990 * t2051;
t2019 = sin(qJ(3,1));
t2006 = 0.1e1 / t2019;
t2100 = t1947 * t2006;
t1948 = -t1981 * t1990 + t1993 * t2051;
t2099 = t1948 * t2006;
t2014 = sin(qJ(2,3));
t1949 = t1979 * t2023 + t2014 * t2053;
t2098 = t1949 * t2004;
t1950 = -t1979 * t2014 + t2023 * t2053;
t2097 = t1950 * t2004;
t2017 = sin(qJ(2,2));
t1951 = t1980 * t2026 + t2017 * t2052;
t2096 = t1951 * t2005;
t1952 = -t1980 * t2017 + t2026 * t2052;
t2095 = t1952 * t2005;
t2020 = sin(qJ(2,1));
t1953 = t1981 * t2029 + t2020 * t2051;
t2094 = t1953 * t2006;
t1954 = -t1981 * t2020 + t2029 * t2051;
t2093 = t1954 * t2006;
t2091 = t1973 * t2024;
t2084 = t1994 * t2015;
t2083 = t1995 * t2018;
t2082 = t1996 * t2021;
t2081 = t1997 * t2015;
t2080 = t1998 * t2018;
t2079 = t1999 * t2021;
t2078 = t2013 * t2014;
t2077 = t2013 * t2023;
t2076 = t2016 * t2017;
t2075 = t2016 * t2026;
t2074 = t2019 * t2020;
t2073 = t2019 * t2029;
t2071 = t1961 * t2091;
t2069 = t1959 * t2089;
t2067 = t1960 * t2087;
t2022 = cos(qJ(3,3));
t1985 = pkin(3) * t2022 + pkin(2);
t1970 = -pkin(3) * t2078 + t1985 * t2023;
t2066 = t1970 * t2004 * t2024;
t2025 = cos(qJ(3,2));
t1986 = pkin(3) * t2025 + pkin(2);
t1971 = -pkin(3) * t2076 + t1986 * t2026;
t2065 = t1971 * t2005 * t2027;
t2028 = cos(qJ(3,1));
t1987 = pkin(3) * t2028 + pkin(2);
t1972 = -pkin(3) * t2074 + t1987 * t2029;
t2064 = t1972 * t2006 * t2030;
t2063 = t1994 * t2091;
t2062 = t1997 * t2091;
t2056 = (cos(qJ(1,3) - t2007) + cos(qJ(1,3) + t2007)) * t2004 / 0.2e1;
t2055 = (cos(qJ(1,2) - t2008) + cos(qJ(1,2) + t2008)) * t2005 / 0.2e1;
t2054 = (cos(qJ(1,1) - t2009) + cos(qJ(1,1) + t2009)) * t2006 / 0.2e1;
t2050 = t1991 * t2071;
t2049 = t2023 * t2071;
t2046 = t2017 * t2069;
t2045 = t2026 * t2069;
t2042 = t2020 * t2067;
t2041 = t2029 * t2067;
t2040 = t1961 * t2063;
t2039 = t1961 * t2062;
t2038 = -t2014 * t2022 - t2077;
t2037 = -t2022 * t2023 + t2078;
t2036 = -t2017 * t2025 - t2075;
t2035 = -t2025 * t2026 + t2076;
t2034 = -t2020 * t2028 - t2073;
t2033 = -t2028 * t2029 + t2074;
t2032 = 0.1e1 / pkin(2);
t2031 = 0.1e1 / pkin(3);
t1969 = pkin(3) * t2073 + t2020 * t1987;
t1968 = pkin(3) * t2075 + t2017 * t1986;
t1967 = pkin(3) * t2077 + t2014 * t1985;
t1942 = t1969 * t1999 + t1972 * t2082;
t1941 = t1968 * t1998 + t1971 * t2083;
t1940 = t1967 * t1997 + t1970 * t2084;
t1939 = t1996 * t1969 - t1972 * t2079;
t1938 = t1995 * t1968 - t1971 * t2080;
t1937 = t1994 * t1967 - t1970 * t2081;
t1936 = t2034 * t1996 - t2033 * t2079;
t1935 = t2034 * t1999 + t2033 * t2082;
t1934 = t2036 * t1995 - t2035 * t2080;
t1933 = t2036 * t1998 + t2035 * t2083;
t1932 = t2038 * t1994 - t2037 * t2081;
t1931 = t2038 * t1997 + t2037 * t2084;
t1 = [0, -t2039 - t2043 - t2047, t2051 * t2058 + t2052 * t2060 + t2053 * t2062, 0, 0, 0, 0, 0, -t1997 * t2049 - t1998 * t2045 - t1999 * t2041 + (t1932 * t2098 + t1934 * t2096 + t1936 * t2094) * t2032, t1998 * t2046 + t1999 * t2042 + t2014 * t2039 + (t1932 * t2097 + t1934 * t2095 + t1936 * t2093) * t2032, 0, 0, 0, 0, 0, -t1997 * t2050 - t1992 * t2047 - t1993 * t2043 + (t1932 * t2104 + t1934 * t2102 + t1936 * t2100 + (t1937 * t2104 + t1938 * t2102 + t1939 * t2100) * t2031) * t2032, t1989 * t2047 + t1990 * t2043 + t1988 * t2039 + (t1932 * t2103 + t1934 * t2101 + t1936 * t2099 + (t1937 * t2103 + t1938 * t2101 + t1939 * t2099) * t2031) * t2032, -g(1); 0, t2040 + t2044 + t2048, -t2051 * t2059 - t2052 * t2061 - t2053 * t2063, 0, 0, 0, 0, 0, t1994 * t2049 + t1995 * t2045 + t1996 * t2041 + (t1931 * t2098 + t1933 * t2096 + t1935 * t2094) * t2032, -t1995 * t2046 - t1996 * t2042 - t2014 * t2040 + (t1931 * t2097 + t1933 * t2095 + t1935 * t2093) * t2032, 0, 0, 0, 0, 0, t1994 * t2050 + t1992 * t2048 + t1993 * t2044 + (t1931 * t2104 + t1933 * t2102 + t1935 * t2100 + (t1940 * t2104 + t1941 * t2102 + t1942 * t2100) * t2031) * t2032, -t1989 * t2048 - t1990 * t2044 - t1988 * t2040 + (t1931 * t2103 + t1933 * t2101 + t1935 * t2099 + (t1940 * t2103 + t1941 * t2101 + t1942 * t2099) * t2031) * t2032, -g(2); 0, t2068 + t2107 + t2070, -t2051 * t2088 - t2052 * t2090 - t2053 * t2092, 0, 0, 0, 0, 0, t2023 * t2107 + t2026 * t2070 + t2029 * t2068 + (t1949 * t2056 + t1951 * t2055 + t1953 * t2054) * t2032, -t2017 * t2070 - t2020 * t2068 - t2014 * t2107 + (t1950 * t2056 + t1952 * t2055 + t1954 * t2054) * t2032, 0, 0, 0, 0, 0, t1991 * t2107 + t1992 * t2070 + t1993 * t2068 + (t1947 * t2054 + t1945 * t2055 + t1943 * t2056 + (-t1943 * t2066 - t1945 * t2065 - t1947 * t2064) * t2031) * t2032, -t1989 * t2070 - t1990 * t2068 - t1988 * t2107 + (t1948 * t2054 + t1946 * t2055 + t1944 * t2056 + (-t1944 * t2066 - t1946 * t2065 - t1948 * t2064) * t2031) * t2032, -g(3);];
tau_reg  = t1;
