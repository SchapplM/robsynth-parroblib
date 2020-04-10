% Calculate Gravitation load for parallel robot
% P3RRRRR2G2P2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2P2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:53
% EndTime: 2020-03-09 21:08:54
% DurationCPUTime: 0.76s
% Computational Cost: add. (630->156), mult. (951->237), div. (69->14), fcn. (603->48), ass. (0->125)
t1021 = cos(qJ(3,1));
t1085 = t1021 ^ 2;
t1018 = cos(qJ(3,2));
t1084 = t1018 ^ 2;
t1015 = cos(qJ(3,3));
t1083 = t1015 ^ 2;
t1012 = sin(qJ(3,1));
t1030 = t1021 * mrSges(3,1) - t1012 * mrSges(3,2);
t1009 = sin(qJ(3,2));
t1031 = t1018 * mrSges(3,1) - t1009 * mrSges(3,2);
t1006 = sin(qJ(3,3));
t1032 = t1015 * mrSges(3,1) - t1006 * mrSges(3,2);
t1082 = -2 * pkin(1);
t1081 = -2 * mrSges(2,1);
t1003 = legFrame(3,2);
t981 = sin(t1003);
t984 = cos(t1003);
t953 = t984 * g(1) - t981 * g(2);
t1080 = mrSges(3,1) * t953;
t1004 = legFrame(2,2);
t982 = sin(t1004);
t985 = cos(t1004);
t954 = t985 * g(1) - t982 * g(2);
t1079 = mrSges(3,1) * t954;
t1005 = legFrame(1,2);
t983 = sin(t1005);
t986 = cos(t1005);
t955 = t986 * g(1) - t983 * g(2);
t1078 = mrSges(3,1) * t955;
t1077 = mrSges(3,2) * t953;
t1076 = mrSges(3,2) * t954;
t1075 = mrSges(3,2) * t955;
t1008 = sin(qJ(1,3));
t1074 = pkin(1) * t1008;
t1011 = sin(qJ(1,2));
t1073 = pkin(1) * t1011;
t1014 = sin(qJ(1,1));
t1072 = pkin(1) * t1014;
t1002 = mrSges(3,3) - mrSges(2,2);
t980 = (t1002 * g(3));
t1017 = cos(qJ(1,3));
t1024 = mrSges(3,2) * g(3);
t1025 = mrSges(1,2) * g(3);
t1026 = mrSges(3,1) * g(3);
t1027 = mrSges(2,1) * g(3);
t967 = (m(2) + m(3)) * pkin(1) + mrSges(1,1);
t959 = t967 * g(3);
t960 = 2 * t980;
t1054 = qJ(2,3) + qJ(3,3);
t974 = qJ(1,3) + t1054;
t961 = cos(t974);
t1055 = qJ(2,3) - qJ(3,3);
t975 = qJ(1,3) + t1055;
t962 = cos(t975);
t999 = qJ(1,3) + qJ(2,3);
t968 = sin(t999);
t971 = cos(t999);
t1007 = sin(qJ(2,3));
t987 = 0.1e1 / t1007;
t1065 = ((-t1024 - t1080) * t962 / 0.2e1 + (t1026 - t1077) * sin(t975) / 0.2e1 + (t1024 - t1080) * t961 / 0.2e1 + (t1026 + t1077) * sin(t974) / 0.2e1 + (t953 * t1081 - t960) * t971 / 0.2e1 + t1008 * (t953 * mrSges(1,2) + t959) + (-t953 * t1002 + t1027) * t968 + (-t953 * t967 + t1025) * t1017) * t987;
t1020 = cos(qJ(1,2));
t1056 = qJ(2,2) + qJ(3,2);
t976 = qJ(1,2) + t1056;
t963 = cos(t976);
t1057 = qJ(2,2) - qJ(3,2);
t977 = qJ(1,2) + t1057;
t964 = cos(t977);
t1000 = qJ(1,2) + qJ(2,2);
t969 = sin(t1000);
t972 = cos(t1000);
t1010 = sin(qJ(2,2));
t988 = 0.1e1 / t1010;
t1064 = ((-t1024 - t1079) * t964 / 0.2e1 + (t1026 - t1076) * sin(t977) / 0.2e1 + (t1024 - t1079) * t963 / 0.2e1 + (t1026 + t1076) * sin(t976) / 0.2e1 + (t954 * t1081 - t960) * t972 / 0.2e1 + t1011 * (t954 * mrSges(1,2) + t959) + (-t954 * t1002 + t1027) * t969 + (-t954 * t967 + t1025) * t1020) * t988;
t1023 = cos(qJ(1,1));
t1058 = qJ(2,1) + qJ(3,1);
t978 = qJ(1,1) + t1058;
t965 = cos(t978);
t1059 = qJ(2,1) - qJ(3,1);
t979 = qJ(1,1) + t1059;
t966 = cos(t979);
t1001 = qJ(1,1) + qJ(2,1);
t970 = sin(t1001);
t973 = cos(t1001);
t1013 = sin(qJ(2,1));
t989 = 0.1e1 / t1013;
t1063 = ((-t1024 - t1078) * t966 / 0.2e1 + (t1026 - t1075) * sin(t979) / 0.2e1 + (t1024 - t1078) * t965 / 0.2e1 + (t1026 + t1075) * sin(t978) / 0.2e1 + (t955 * t1081 - t960) * t973 / 0.2e1 + t1014 * (t955 * mrSges(1,2) + t959) + (-t955 * t1002 + t1027) * t970 + (-t955 * t967 + t1025) * t1023) * t989;
t991 = 0.1e1 / t1015;
t1062 = t991 * (-(t981 * g(1) + t984 * g(2)) * t1032 + (g(3) * t971 + t953 * t968) * (mrSges(3,1) * t1006 + mrSges(3,2) * t1015));
t994 = 0.1e1 / t1018;
t1061 = t994 * (-(t982 * g(1) + t985 * g(2)) * t1031 + (g(3) * t972 + t954 * t969) * (mrSges(3,1) * t1009 + mrSges(3,2) * t1018));
t997 = 0.1e1 / t1021;
t1060 = t997 * (-(t983 * g(1) + t986 * g(2)) * t1030 + (g(3) * t973 + t955 * t970) * (mrSges(3,1) * t1012 + mrSges(3,2) * t1021));
t1016 = cos(qJ(2,3));
t950 = t1017 * t1007 + t1008 * t1016;
t1053 = t1015 * t950;
t1019 = cos(qJ(2,2));
t951 = t1020 * t1010 + t1011 * t1019;
t1052 = t1018 * t951;
t1022 = cos(qJ(2,1));
t952 = t1023 * t1013 + t1014 * t1022;
t1051 = t1021 * t952;
t1050 = t981 * t1006;
t1049 = t982 * t1009;
t1048 = t983 * t1012;
t1047 = t984 * t1006;
t1046 = t985 * t1009;
t1045 = t986 * t1012;
t1044 = t950 * t1083 * pkin(2);
t1043 = t951 * t1084 * pkin(2);
t1042 = t952 * t1085 * pkin(2);
t1041 = t1016 * t1006 * pkin(1);
t1040 = t1019 * t1009 * pkin(1);
t1039 = t1022 * t1012 * pkin(1);
t1038 = t991 * t1065;
t1037 = t994 * t1064;
t1036 = t997 * t1063;
t944 = -t980 * t971 + (t1032 * g(3) + t1027) * t968 + ((-mrSges(2,1) - t1032) * t971 - t1002 * t968) * t953;
t1035 = t944 * t987 / t1083;
t945 = -t980 * t972 + (t1031 * g(3) + t1027) * t969 + ((-mrSges(2,1) - t1031) * t972 - t1002 * t969) * t954;
t1034 = t945 * t988 / t1084;
t946 = -t980 * t973 + (t1030 * g(3) + t1027) * t970 + ((-mrSges(2,1) - t1030) * t973 - t1002 * t970) * t955;
t1033 = t946 * t989 / t1085;
t1029 = 1 / pkin(1);
t1028 = 0.1e1 / pkin(2);
t1 = [-g(1) * m(4) + ((t986 * t1051 + t1048) * t1036 + (t985 * t1052 + t1049) * t1037 + (t984 * t1053 + t1050) * t1038) * t1029 + (t981 * t1062 + t982 * t1061 + t983 * t1060 + ((-t986 * t1042 + (-pkin(2) * t1048 - t986 * t1072) * t1021 - t983 * t1039) * t1033 + (-t985 * t1043 + (-pkin(2) * t1049 - t985 * t1073) * t1018 - t982 * t1040) * t1034 + (-t984 * t1044 + (-pkin(2) * t1050 - t984 * t1074) * t1015 - t981 * t1041) * t1035) * t1029) * t1028; -g(2) * m(4) + ((-t983 * t1051 + t1045) * t1036 + (-t982 * t1052 + t1046) * t1037 + (-t981 * t1053 + t1047) * t1038) * t1029 + (t984 * t1062 + t985 * t1061 + t986 * t1060 + ((t983 * t1042 + (-pkin(2) * t1045 + t983 * t1072) * t1021 - t986 * t1039) * t1033 + (t982 * t1043 + (-pkin(2) * t1046 + t982 * t1073) * t1018 - t985 * t1040) * t1034 + (t981 * t1044 + (-pkin(2) * t1047 + t981 * t1074) * t1015 - t984 * t1041) * t1035) * t1029) * t1028; -g(3) * m(4) + (t971 * t1065 + t972 * t1064 + t973 * t1063 + ((t1023 * t1082 + (-t965 - t966) * pkin(2)) / (sin(t1058) + sin(t1059)) * t946 + (t1020 * t1082 + (-t963 - t964) * pkin(2)) / (sin(t1056) + sin(t1057)) * t945 + (t1017 * t1082 + (-t961 - t962) * pkin(2)) / (sin(t1054) + sin(t1055)) * t944) * t1028) * t1029;];
taugX  = t1;
