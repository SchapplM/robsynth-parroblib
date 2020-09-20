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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
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

function taugX = P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:47:12
% EndTime: 2020-08-06 18:47:13
% DurationCPUTime: 0.50s
% Computational Cost: add. (585->125), mult. (687->182), div. (30->7), fcn. (492->32), ass. (0->90)
t1031 = 2 * pkin(1);
t1006 = (m(2) + m(3));
t1030 = mrSges(3,1) * g(3);
t982 = pkin(7) + qJ(3,3);
t967 = cos(t982);
t1029 = pkin(3) * t967;
t983 = pkin(7) + qJ(3,2);
t968 = cos(t983);
t1028 = pkin(3) * t968;
t984 = pkin(7) + qJ(3,1);
t969 = cos(t984);
t1027 = pkin(3) * t969;
t992 = pkin(5) + qJ(2,1);
t991 = pkin(5) + qJ(2,2);
t990 = pkin(5) + qJ(2,3);
t1000 = cos(qJ(1,3));
t987 = legFrame(3,3);
t970 = sin(t987);
t973 = cos(t987);
t950 = -t970 * g(1) + t973 * g(2);
t953 = t973 * g(1) + t970 * g(2);
t979 = -pkin(6) - t990;
t976 = 0.1e1 / t979;
t994 = sin(qJ(1,3));
t1025 = (t1000 * t950 - t953 * t994) * t976;
t1002 = cos(qJ(1,2));
t988 = legFrame(2,3);
t971 = sin(t988);
t974 = cos(t988);
t951 = -t971 * g(1) + t974 * g(2);
t954 = t974 * g(1) + t971 * g(2);
t980 = -pkin(6) - t991;
t977 = 0.1e1 / t980;
t996 = sin(qJ(1,2));
t1024 = (t1002 * t951 - t954 * t996) * t977;
t1004 = cos(qJ(1,1));
t989 = legFrame(1,3);
t972 = sin(t989);
t975 = cos(t989);
t952 = -t972 * g(1) + t975 * g(2);
t955 = t975 * g(1) + t972 * g(2);
t981 = -pkin(6) - t992;
t978 = 0.1e1 / t981;
t998 = sin(qJ(1,1));
t1023 = (t1004 * t952 - t955 * t998) * t978;
t1019 = mrSges(2,3) + mrSges(3,3) - mrSges(1,2);
t1009 = m(2) * qJ(2,3) + t990 * m(3) + t1019;
t1018 = -m(3) * pkin(2) - mrSges(2,1);
t959 = pkin(1) * t1006 + mrSges(1,1);
t985 = sin(pkin(7));
t986 = cos(pkin(7));
t993 = sin(qJ(3,3));
t999 = cos(qJ(3,3));
t1014 = (-mrSges(3,1) * t999 + mrSges(3,2) * t993 + t1018) * t986 + (t993 * mrSges(3,1) + mrSges(3,2) * t999 + mrSges(2,2)) * t985 - t959;
t1022 = t976 * ((-t1009 * t1000 - t1014 * t994) * t953 + (t1014 * t1000 - t1009 * t994) * t950);
t1010 = m(2) * qJ(2,2) + t991 * m(3) + t1019;
t1001 = cos(qJ(3,2));
t995 = sin(qJ(3,2));
t1013 = (-mrSges(3,1) * t1001 + mrSges(3,2) * t995 + t1018) * t986 + (t995 * mrSges(3,1) + mrSges(3,2) * t1001 + mrSges(2,2)) * t985 - t959;
t1021 = t977 * ((-t1010 * t1002 - t1013 * t996) * t954 + (t1013 * t1002 - t1010 * t996) * t951);
t1011 = m(2) * qJ(2,1) + t992 * m(3) + t1019;
t1003 = cos(qJ(3,1));
t997 = sin(qJ(3,1));
t1012 = (-mrSges(3,1) * t1003 + mrSges(3,2) * t997 + t1018) * t986 + (t997 * mrSges(3,1) + mrSges(3,2) * t1003 + mrSges(2,2)) * t985 - t959;
t1020 = t978 * ((-t1011 * t1004 - t1012 * t998) * t955 + (t1012 * t1004 - t1011 * t998) * t952);
t1017 = t953 * t1000 + t950 * t994;
t1016 = t954 * t1002 + t951 * t996;
t1015 = t955 * t1004 + t952 * t998;
t1007 = 0.2e1 * pkin(7);
t1005 = mrSges(3,2) * g(3);
t966 = sin(t984);
t965 = sin(t983);
t964 = sin(t982);
t963 = 0.1e1 / t969;
t962 = 0.1e1 / t968;
t961 = 0.1e1 / t967;
t960 = t986 * pkin(2) + pkin(1);
t946 = t972 * t1004 + t975 * t998;
t945 = t971 * t1002 + t974 * t996;
t944 = t970 * t1000 + t973 * t994;
t943 = t975 * t1004 - t972 * t998;
t942 = t974 * t1002 - t971 * t996;
t941 = t973 * t1000 - t970 * t994;
t940 = t960 * t1004 - t998 * t981;
t939 = t960 * t1002 - t996 * t980;
t938 = t960 * t1000 - t994 * t979;
t937 = t1004 * t981 + t998 * t960;
t936 = t1002 * t980 + t996 * t960;
t935 = t1000 * t979 + t994 * t960;
t1 = [-t941 * t1022 - t942 * t1021 - t943 * t1020 - g(1) * m(4) + (-(t943 * t1027 - t937 * t972 + t940 * t975) * t1023 - (t942 * t1028 - t936 * t971 + t939 * t974) * t1024 - (t941 * t1029 - t935 * t970 + t938 * t973) * t1025) * t1006; -t944 * t1022 - t945 * t1021 - t946 * t1020 - g(2) * m(4) + (-(t946 * t1027 + t937 * t975 + t940 * t972) * t1023 - (t945 * t1028 + t936 * t974 + t939 * t971) * t1024 - (t944 * t1029 + t935 * t973 + t938 * t970) * t1025) * t1006; -t966 * t963 * t1020 - t965 * t962 * t1021 - t964 * t961 * t1022 - g(3) * m(4) - ((pkin(3) * sin(0.2e1 * t984) + t966 * t1031 + (sin(t1007 + qJ(3,1)) + t997) * pkin(2)) * t963 * t1023 + (pkin(3) * sin(0.2e1 * t983) + t965 * t1031 + (sin(t1007 + qJ(3,2)) + t995) * pkin(2)) * t962 * t1024 + (pkin(3) * sin(0.2e1 * t982) + t964 * t1031 + (sin(t1007 + qJ(3,3)) + t993) * pkin(2)) * t961 * t1025) * t1006 / 0.2e1 + (t963 * ((t1015 * mrSges(3,2) - t1030) * t969 + t966 * (t1015 * mrSges(3,1) + t1005)) + t962 * ((t1016 * mrSges(3,2) - t1030) * t968 + t965 * (t1016 * mrSges(3,1) + t1005)) + t961 * ((t1017 * mrSges(3,2) - t1030) * t967 + t964 * (t1017 * mrSges(3,1) + t1005))) / pkin(3);];
taugX  = t1;
