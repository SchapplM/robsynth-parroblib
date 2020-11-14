% Calculate Gravitation load for parallel robot
% P3RPRRR9V1G1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:46:54
% EndTime: 2020-08-06 18:46:55
% DurationCPUTime: 0.72s
% Computational Cost: add. (621->132), mult. (924->187), div. (30->7), fcn. (564->32), ass. (0->90)
t967 = m(2) + m(3);
t996 = 2 * pkin(1);
t995 = m(2) * rSges(2,1);
t994 = rSges(3,1) * g(3);
t993 = -t967 / 0.2e1;
t992 = m(3) / pkin(3);
t944 = pkin(7) + qJ(3,3);
t926 = cos(t944);
t991 = pkin(3) * t926;
t945 = pkin(7) + qJ(3,2);
t927 = cos(t945);
t990 = pkin(3) * t927;
t946 = pkin(7) + qJ(3,1);
t928 = cos(t946);
t989 = pkin(3) * t928;
t988 = pkin(5) + qJ(2,1);
t987 = pkin(5) + qJ(2,2);
t986 = pkin(5) + qJ(2,3);
t949 = legFrame(3,3);
t929 = sin(t949);
t932 = cos(t949);
t909 = -t929 * g(1) + t932 * g(2);
t912 = t932 * g(1) + t929 * g(2);
t941 = -pkin(6) - t986;
t935 = 0.1e1 / t941;
t953 = sin(qJ(1,3));
t959 = cos(qJ(1,3));
t985 = (t909 * t959 - t912 * t953) * t935;
t950 = legFrame(2,3);
t930 = sin(t950);
t933 = cos(t950);
t910 = -t930 * g(1) + t933 * g(2);
t913 = t933 * g(1) + t930 * g(2);
t942 = -pkin(6) - t987;
t936 = 0.1e1 / t942;
t955 = sin(qJ(1,2));
t961 = cos(qJ(1,2));
t984 = (t910 * t961 - t913 * t955) * t936;
t951 = legFrame(1,3);
t931 = sin(t951);
t934 = cos(t951);
t911 = -t931 * g(1) + t934 * g(2);
t914 = t934 * g(1) + t931 * g(2);
t943 = -pkin(6) - t988;
t937 = 0.1e1 / t943;
t957 = sin(qJ(1,1));
t963 = cos(qJ(1,1));
t983 = (t911 * t963 - t914 * t957) * t937;
t947 = sin(pkin(7));
t948 = cos(pkin(7));
t952 = sin(qJ(3,3));
t958 = cos(qJ(3,3));
t965 = rSges(2,2) * m(2);
t979 = m(1) * rSges(1,1) + pkin(1) * t967;
t972 = -(-t995 + (-rSges(3,1) * t958 + rSges(3,2) * t952 - pkin(2)) * m(3)) * t948 - (t965 + (rSges(3,1) * t952 + rSges(3,2) * t958) * m(3)) * t947 + t979;
t966 = m(1) * rSges(1,2);
t978 = -(rSges(3,3) + t986) * m(3) + (-rSges(2,3) - qJ(2,3)) * m(2) + t966;
t982 = t935 * ((t972 * t953 + t978 * t959) * t912 + (t978 * t953 - t972 * t959) * t909);
t954 = sin(qJ(3,2));
t960 = cos(qJ(3,2));
t971 = -(-t995 + (-rSges(3,1) * t960 + rSges(3,2) * t954 - pkin(2)) * m(3)) * t948 - (t965 + (rSges(3,1) * t954 + rSges(3,2) * t960) * m(3)) * t947 + t979;
t977 = -(rSges(3,3) + t987) * m(3) + (-rSges(2,3) - qJ(2,2)) * m(2) + t966;
t981 = t936 * ((t971 * t955 + t977 * t961) * t913 + (t977 * t955 - t971 * t961) * t910);
t956 = sin(qJ(3,1));
t962 = cos(qJ(3,1));
t970 = -(-t995 + (-rSges(3,1) * t962 + rSges(3,2) * t956 - pkin(2)) * m(3)) * t948 - (t965 + (rSges(3,1) * t956 + rSges(3,2) * t962) * m(3)) * t947 + t979;
t976 = -(rSges(3,3) + t988) * m(3) + (-rSges(2,3) - qJ(2,1)) * m(2) + t966;
t980 = t937 * ((t970 * t957 + t976 * t963) * t914 + (t976 * t957 - t970 * t963) * t911);
t975 = t909 * t953 + t912 * t959;
t974 = t910 * t955 + t913 * t961;
t973 = t911 * t957 + t914 * t963;
t968 = 0.2e1 * pkin(7);
t964 = rSges(3,2) * g(3);
t925 = sin(t946);
t924 = sin(t945);
t923 = sin(t944);
t918 = t948 * pkin(2) + pkin(1);
t905 = t931 * t963 + t934 * t957;
t904 = t930 * t961 + t933 * t955;
t903 = t929 * t959 + t932 * t953;
t902 = -t931 * t957 + t934 * t963;
t901 = -t930 * t955 + t933 * t961;
t900 = -t929 * t953 + t932 * t959;
t896 = t918 * t963 - t957 * t943;
t895 = t918 * t961 - t955 * t942;
t894 = t918 * t959 - t953 * t941;
t893 = t957 * t918 + t963 * t943;
t892 = t955 * t918 + t961 * t942;
t891 = t953 * t918 + t959 * t941;
t1 = [-t900 * t982 - t901 * t981 - t902 * t980 - m(4) * g(1) + (-(-t893 * t931 + t896 * t934 + t902 * t989) * t983 - (-t892 * t930 + t895 * t933 + t901 * t990) * t984 - (-t891 * t929 + t894 * t932 + t900 * t991) * t985) * t967; -t903 * t982 - t904 * t981 - t905 * t980 - m(4) * g(2) + (-(t893 * t934 + t896 * t931 + t905 * t989) * t983 - (t892 * t933 + t895 * t930 + t904 * t990) * t984 - (t891 * t932 + t894 * t929 + t903 * t991) * t985) * t967; -m(4) * g(3) + (-t925 * t980 + (pkin(3) * sin(0.2e1 * t946) + t925 * t996 + (sin(t968 + qJ(3,1)) + t956) * pkin(2)) * t983 * t993 + ((t973 * rSges(3,2) - t994) * t928 + t925 * (t973 * rSges(3,1) + t964)) * t992) / t928 + (-t924 * t981 + (pkin(3) * sin(0.2e1 * t945) + t924 * t996 + (sin(t968 + qJ(3,2)) + t954) * pkin(2)) * t984 * t993 + ((t974 * rSges(3,2) - t994) * t927 + t924 * (t974 * rSges(3,1) + t964)) * t992) / t927 + (-t923 * t982 + (pkin(3) * sin(0.2e1 * t944) + t923 * t996 + (sin(t968 + qJ(3,3)) + t952) * pkin(2)) * t985 * t993 + ((t975 * rSges(3,2) - t994) * t926 + t923 * (t975 * rSges(3,1) + t964)) * t992) / t926;];
taugX  = t1;
