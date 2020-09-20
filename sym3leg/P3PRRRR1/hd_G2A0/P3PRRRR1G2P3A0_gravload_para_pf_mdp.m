% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR1G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G2P3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G2P3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G2P3A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:46
% EndTime: 2020-03-09 21:16:47
% DurationCPUTime: 0.67s
% Computational Cost: add. (192->66), mult. (444->136), div. (144->10), fcn. (504->18), ass. (0->62)
t1942 = legFrame(3,2);
t1921 = sin(t1942);
t1924 = cos(t1942);
t1915 = g(1) * t1921 + g(2) * t1924;
t1946 = sin(qJ(2,3));
t1952 = cos(qJ(2,3));
t1906 = -g(3) * t1946 + t1915 * t1952;
t1945 = sin(qJ(3,3));
t2008 = t1906 * t1945;
t1943 = legFrame(2,2);
t1922 = sin(t1943);
t1925 = cos(t1943);
t1916 = g(1) * t1922 + g(2) * t1925;
t1948 = sin(qJ(2,2));
t1954 = cos(qJ(2,2));
t1907 = -g(3) * t1948 + t1916 * t1954;
t1947 = sin(qJ(3,2));
t2007 = t1907 * t1947;
t1944 = legFrame(1,2);
t1923 = sin(t1944);
t1926 = cos(t1944);
t1917 = g(1) * t1923 + g(2) * t1926;
t1950 = sin(qJ(2,1));
t1956 = cos(qJ(2,1));
t1908 = -g(3) * t1950 + t1917 * t1956;
t1949 = sin(qJ(3,1));
t2006 = t1908 * t1949;
t1920 = g(1) * t1926 - g(2) * t1923;
t1955 = cos(qJ(3,1));
t1940 = 0.1e1 / t1955;
t1967 = g(3) * t1956 + t1917 * t1950;
t1973 = t1940 * t2006;
t1935 = 0.1e1 / t1950;
t1982 = t1935 * t1956;
t2005 = ((MDP(11) * t2006 - MDP(3) * t1908 + MDP(4) * t1967) * t1949 / t1955 ^ 2 - MDP(10) * t1973) * t1982 - (MDP(10) * (t1920 * t1955 + t1949 * t1967) + (-t1920 * t1949 + t1955 * t1967) * MDP(11)) * t1940;
t1919 = g(1) * t1925 - g(2) * t1922;
t1953 = cos(qJ(3,2));
t1938 = 0.1e1 / t1953;
t1968 = g(3) * t1954 + t1916 * t1948;
t1974 = t1938 * t2007;
t1933 = 0.1e1 / t1948;
t1984 = t1933 * t1954;
t2004 = ((MDP(11) * t2007 - MDP(3) * t1907 + MDP(4) * t1968) * t1947 / t1953 ^ 2 - MDP(10) * t1974) * t1984 - (MDP(10) * (t1919 * t1953 + t1947 * t1968) + (-t1919 * t1947 + t1953 * t1968) * MDP(11)) * t1938;
t1918 = g(1) * t1924 - g(2) * t1921;
t1951 = cos(qJ(3,3));
t1936 = 0.1e1 / t1951;
t1969 = g(3) * t1952 + t1915 * t1946;
t1975 = t1936 * t2008;
t1931 = 0.1e1 / t1946;
t1986 = t1931 * t1952;
t2003 = ((MDP(11) * t2008 - MDP(3) * t1906 + MDP(4) * t1969) * t1945 / t1951 ^ 2 - MDP(10) * t1975) * t1986 - (MDP(10) * (t1918 * t1951 + t1945 * t1969) + (-t1918 * t1945 + t1951 * t1969) * MDP(11)) * t1936;
t1987 = t1931 * t1936;
t1985 = t1933 * t1938;
t1983 = t1935 * t1940;
t1981 = t1946 * t1951;
t1980 = t1948 * t1953;
t1979 = t1950 * t1955;
t1972 = t1915 * t1987;
t1971 = t1916 * t1985;
t1970 = t1917 * t1983;
t1957 = 0.1e1 / pkin(2);
t1 = [(-(t1923 * t1979 - t1926 * t1949) * t1970 - (t1922 * t1980 - t1925 * t1947) * t1971 - (t1921 * t1981 - t1924 * t1945) * t1972) * MDP(1) - g(1) * MDP(12) + (t2003 * t1924 + t2004 * t1925 + t2005 * t1926) * t1957; (-(t1923 * t1949 + t1926 * t1979) * t1970 - (t1922 * t1947 + t1925 * t1980) * t1971 - (t1921 * t1945 + t1924 * t1981) * t1972) * MDP(1) - g(2) * MDP(12) + (-t2003 * t1921 - t2004 * t1922 - t2005 * t1923) * t1957; (-t1915 * t1986 - t1916 * t1984 - t1917 * t1982) * MDP(1) - g(3) * MDP(12) + ((t1906 * t1987 + t1907 * t1985 + t1908 * t1983) * MDP(3) + (-t1967 * t1983 - t1968 * t1985 - t1969 * t1987) * MDP(4) + (t1906 * t1931 + t1907 * t1933 + t1908 * t1935) * MDP(10) + (-t1931 * t1975 - t1933 * t1974 - t1935 * t1973) * MDP(11)) * t1957;];
taugX  = t1;
