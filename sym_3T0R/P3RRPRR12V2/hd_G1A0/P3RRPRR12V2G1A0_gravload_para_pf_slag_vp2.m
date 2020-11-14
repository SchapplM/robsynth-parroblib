% Calculate Gravitation load for parallel robot
% P3RRPRR12V2G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:15:13
% EndTime: 2020-08-06 19:15:14
% DurationCPUTime: 0.72s
% Computational Cost: add. (717->154), mult. (909->247), div. (36->6), fcn. (636->18), ass. (0->114)
t1118 = m(2) + m(3);
t1076 = pkin(2) + pkin(3);
t1117 = mrSges(3,3) - mrSges(2,2);
t1044 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t1116 = g(3) * t1044;
t1063 = sin(qJ(2,3));
t1115 = t1063 * qJ(3,3);
t1065 = sin(qJ(2,2));
t1114 = t1065 * qJ(3,2);
t1067 = sin(qJ(2,1));
t1113 = t1067 * qJ(3,1);
t1041 = m(3) * qJ(3,3) + t1117;
t1069 = cos(qJ(2,3));
t1077 = 0.1e1 / qJ(3,3);
t1060 = legFrame(3,3);
t1045 = sin(t1060);
t1048 = cos(t1060);
t1015 = -t1045 * g(1) + t1048 * g(2);
t1018 = t1048 * g(1) + t1045 * g(2);
t1064 = sin(qJ(1,3));
t1070 = cos(qJ(1,3));
t1088 = t1015 * t1064 + t1018 * t1070;
t1112 = t1077 * ((-t1088 * t1041 - t1116) * t1069 + (-g(3) * t1041 + t1088 * t1044) * t1063);
t1111 = t1077 * (-t1069 * g(3) + t1088 * t1063);
t1042 = m(3) * qJ(3,2) + t1117;
t1071 = cos(qJ(2,2));
t1078 = 0.1e1 / qJ(3,2);
t1061 = legFrame(2,3);
t1046 = sin(t1061);
t1049 = cos(t1061);
t1016 = -t1046 * g(1) + t1049 * g(2);
t1019 = t1049 * g(1) + t1046 * g(2);
t1066 = sin(qJ(1,2));
t1072 = cos(qJ(1,2));
t1087 = t1016 * t1066 + t1019 * t1072;
t1110 = t1078 * ((-t1087 * t1042 - t1116) * t1071 + (-g(3) * t1042 + t1087 * t1044) * t1065);
t1109 = t1078 * (-t1071 * g(3) + t1087 * t1065);
t1043 = qJ(3,1) * m(3) + t1117;
t1073 = cos(qJ(2,1));
t1079 = 0.1e1 / qJ(3,1);
t1062 = legFrame(1,3);
t1047 = sin(t1062);
t1050 = cos(t1062);
t1017 = -t1047 * g(1) + t1050 * g(2);
t1020 = t1050 * g(1) + t1047 * g(2);
t1068 = sin(qJ(1,1));
t1074 = cos(qJ(1,1));
t1086 = t1017 * t1068 + t1020 * t1074;
t1108 = t1079 * ((-t1086 * t1043 - t1116) * t1073 + (-g(3) * t1043 + t1086 * t1044) * t1067);
t1107 = t1079 * (-t1073 * g(3) + t1086 * t1067);
t1075 = pkin(5) - pkin(6);
t1035 = t1064 * t1075;
t1036 = t1066 * t1075;
t1037 = t1068 * t1075;
t1038 = t1075 * t1070;
t1039 = t1075 * t1072;
t1040 = t1075 * t1074;
t1106 = t1076 * t1069;
t1105 = t1076 * t1071;
t1104 = t1076 * t1073;
t1027 = pkin(1) + t1115;
t1021 = 0.1e1 / (t1027 + t1106);
t1103 = t1021 * t1111;
t1029 = pkin(1) + t1114;
t1022 = 0.1e1 / (t1029 + t1105);
t1102 = t1022 * t1109;
t1031 = pkin(1) + t1113;
t1023 = 0.1e1 / (t1031 + t1104);
t1101 = t1023 * t1107;
t1100 = t1069 * t1112;
t1099 = t1071 * t1110;
t1098 = t1073 * t1108;
t1097 = (qJ(3,3) + t1076) * (-qJ(3,3) + t1076) * t1069 ^ 2;
t1096 = (qJ(3,2) + t1076) * (-qJ(3,2) + t1076) * t1071 ^ 2;
t1095 = (qJ(3,1) + t1076) * (-qJ(3,1) + t1076) * t1073 ^ 2;
t1026 = pkin(1) + 0.2e1 * t1115;
t1094 = t1026 * t1070 + t1035;
t1093 = t1027 * t1070 + t1035;
t1028 = pkin(1) + 0.2e1 * t1114;
t1092 = t1028 * t1072 + t1036;
t1091 = t1029 * t1072 + t1036;
t1030 = pkin(1) + 0.2e1 * t1113;
t1090 = t1030 * t1074 + t1037;
t1089 = t1031 * t1074 + t1037;
t1032 = pkin(1) * t1063 + qJ(3,3);
t1085 = t1032 * t1070 + t1063 * t1035;
t1033 = pkin(1) * t1065 + qJ(3,2);
t1084 = t1033 * t1072 + t1065 * t1036;
t1034 = pkin(1) * t1067 + qJ(3,1);
t1083 = t1034 * t1074 + t1067 * t1037;
t1025 = t1118 * pkin(1) + mrSges(1,1);
t1082 = t1041 * t1063 + t1044 * t1069 + t1025;
t1081 = t1042 * t1065 + t1044 * t1071 + t1025;
t1080 = t1043 * t1067 + t1044 * t1073 + t1025;
t1024 = t1118 * pkin(5) - mrSges(1,2) + mrSges(3,2) + mrSges(2,3);
t1014 = t1068 * t1031 - t1040;
t1013 = t1068 * t1030 - t1040;
t1012 = t1066 * t1029 - t1039;
t1011 = t1066 * t1028 - t1039;
t1010 = t1064 * t1027 - t1038;
t1009 = t1064 * t1026 - t1038;
t1008 = t1047 * t1074 + t1050 * t1068;
t1007 = t1047 * t1068 - t1050 * t1074;
t1006 = t1046 * t1072 + t1049 * t1066;
t1005 = t1046 * t1066 - t1049 * t1072;
t1004 = t1045 * t1070 + t1048 * t1064;
t1003 = t1045 * t1064 - t1048 * t1070;
t1002 = t1068 * t1034 - t1067 * t1040;
t1001 = t1066 * t1033 - t1065 * t1039;
t1000 = t1064 * t1032 - t1063 * t1038;
t993 = (-t1024 * t1074 + t1068 * t1080) * t1020 + (-t1068 * t1024 - t1080 * t1074) * t1017;
t992 = (-t1024 * t1072 + t1066 * t1081) * t1019 + (-t1066 * t1024 - t1081 * t1072) * t1016;
t991 = (-t1024 * t1070 + t1064 * t1082) * t1018 + (-t1064 * t1024 - t1082 * t1070) * t1015;
t1 = [-g(1) * m(4) + (-t1008 * t993 - (t1007 * t1104 + t1014 * t1047 - t1089 * t1050) * t1098) * t1023 + (-t1006 * t992 - (t1005 * t1105 + t1012 * t1046 - t1091 * t1049) * t1099) * t1022 + (-t1004 * t991 - (t1003 * t1106 + t1010 * t1045 - t1093 * t1048) * t1100) * t1021 + (-(-t1007 * t1095 - (t1047 * t1013 - t1090 * t1050) * t1104 - (t1047 * t1002 - t1083 * t1050) * qJ(3,1)) * t1101 - (-t1005 * t1096 - (t1046 * t1011 - t1092 * t1049) * t1105 - (t1046 * t1001 - t1084 * t1049) * qJ(3,2)) * t1102 - (-t1003 * t1097 - (t1045 * t1009 - t1094 * t1048) * t1106 - (t1045 * t1000 - t1085 * t1048) * qJ(3,3)) * t1103) * m(3); -g(2) * m(4) + (-t1007 * t993 + (t1008 * t1104 + t1014 * t1050 + t1047 * t1089) * t1098) * t1023 + (-t1005 * t992 + (t1006 * t1105 + t1012 * t1049 + t1046 * t1091) * t1099) * t1022 + (-t1003 * t991 + (t1004 * t1106 + t1010 * t1048 + t1045 * t1093) * t1100) * t1021 + (-(t1008 * t1095 + (t1013 * t1050 + t1090 * t1047) * t1104 + (t1002 * t1050 + t1083 * t1047) * qJ(3,1)) * t1101 - (t1006 * t1096 + (t1011 * t1049 + t1092 * t1046) * t1105 + (t1001 * t1049 + t1084 * t1046) * qJ(3,2)) * t1102 - (t1004 * t1097 + (t1009 * t1048 + t1094 * t1045) * t1106 + (t1000 * t1048 + t1085 * t1045) * qJ(3,3)) * t1103) * m(3); t1063 * t1112 + t1065 * t1110 + t1067 * t1108 - g(3) * m(4) + (-(-t1073 * qJ(3,1) + t1076 * t1067) * t1107 - (-t1071 * qJ(3,2) + t1076 * t1065) * t1109 - (-t1069 * qJ(3,3) + t1076 * t1063) * t1111) * m(3);];
taugX  = t1;
