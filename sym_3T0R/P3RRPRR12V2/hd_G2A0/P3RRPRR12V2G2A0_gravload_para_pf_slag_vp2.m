% Calculate Gravitation load for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:21:06
% EndTime: 2020-08-06 19:21:07
% DurationCPUTime: 0.86s
% Computational Cost: add. (870->178), mult. (1047->312), div. (45->6), fcn. (648->18), ass. (0->134)
t1152 = mrSges(3,3) - mrSges(2,2);
t1046 = qJ(3,1) * m(3) + t1152;
t1047 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t1076 = sin(qJ(2,1));
t1082 = cos(qJ(2,1));
t1157 = t1076 * t1046 + t1047 * t1082;
t1045 = m(3) * qJ(3,2) + t1152;
t1074 = sin(qJ(2,2));
t1080 = cos(qJ(2,2));
t1156 = t1074 * t1045 + t1047 * t1080;
t1044 = m(3) * qJ(3,3) + t1152;
t1072 = sin(qJ(2,3));
t1078 = cos(qJ(2,3));
t1155 = t1072 * t1044 + t1047 * t1078;
t1085 = (pkin(2) + pkin(3));
t1154 = -2 * t1085;
t1153 = m(2) + m(3);
t1069 = legFrame(3,2);
t1051 = sin(t1069);
t1151 = t1051 * qJ(3,3);
t1070 = legFrame(2,2);
t1052 = sin(t1070);
t1150 = t1052 * qJ(3,2);
t1071 = legFrame(1,2);
t1053 = sin(t1071);
t1149 = t1053 * qJ(3,1);
t1054 = cos(t1069);
t1148 = t1054 * qJ(3,3);
t1055 = cos(t1070);
t1147 = t1055 * qJ(3,2);
t1056 = cos(t1071);
t1146 = t1056 * qJ(3,1);
t1073 = sin(qJ(1,3));
t1145 = t1073 * qJ(3,3);
t1075 = sin(qJ(1,2));
t1144 = t1075 * qJ(3,2);
t1077 = sin(qJ(1,1));
t1143 = t1077 * qJ(3,1);
t1028 = t1051 * g(1) + t1054 * g(2);
t1031 = t1054 * g(1) - t1051 * g(2);
t1079 = cos(qJ(1,3));
t1091 = g(3) * t1079 + t1031 * t1073;
t1007 = (-t1028 * t1047 - t1091 * t1044) * t1078 + t1072 * (-t1028 * t1044 + t1091 * t1047);
t1086 = 0.1e1 / qJ(3,3);
t1142 = t1007 * t1086;
t1029 = t1052 * g(1) + t1055 * g(2);
t1032 = t1055 * g(1) - t1052 * g(2);
t1081 = cos(qJ(1,2));
t1090 = g(3) * t1081 + t1032 * t1075;
t1008 = (-t1029 * t1047 - t1090 * t1045) * t1080 + t1074 * (-t1029 * t1045 + t1090 * t1047);
t1087 = 0.1e1 / qJ(3,2);
t1141 = t1008 * t1087;
t1030 = t1053 * g(1) + t1056 * g(2);
t1033 = t1056 * g(1) - t1053 * g(2);
t1083 = cos(qJ(1,1));
t1089 = g(3) * t1083 + t1033 * t1077;
t1009 = (-t1030 * t1047 - t1089 * t1046) * t1082 + t1076 * (-t1030 * t1046 + t1089 * t1047);
t1088 = 0.1e1 / qJ(3,1);
t1140 = t1009 * t1088;
t1084 = pkin(5) - pkin(6);
t1097 = pkin(1) * t1073 - t1084 * t1079;
t1106 = t1072 * t1145;
t1139 = (t1097 + 0.2e1 * t1106) * t1085;
t1096 = pkin(1) * t1075 - t1084 * t1081;
t1105 = t1074 * t1144;
t1138 = (t1096 + 0.2e1 * t1105) * t1085;
t1095 = pkin(1) * t1077 - t1084 * t1083;
t1104 = t1076 * t1143;
t1137 = (t1095 + 0.2e1 * t1104) * t1085;
t1048 = t1072 * qJ(3,3);
t1112 = t1085 * t1078;
t1094 = t1048 + pkin(1) + t1112;
t1025 = 0.1e1 / t1094;
t1136 = t1025 * t1086;
t1049 = t1074 * qJ(3,2);
t1111 = t1085 * t1080;
t1093 = t1049 + pkin(1) + t1111;
t1026 = 0.1e1 / t1093;
t1135 = t1026 * t1087;
t1050 = t1076 * qJ(3,1);
t1110 = t1085 * t1082;
t1092 = t1050 + pkin(1) + t1110;
t1027 = 0.1e1 / t1092;
t1134 = t1027 * t1088;
t1130 = (qJ(3,3) + t1085) * (-qJ(3,3) + t1085);
t1129 = (qJ(3,2) + t1085) * (-qJ(3,2) + t1085);
t1128 = (qJ(3,1) + t1085) * (-qJ(3,1) + t1085);
t1126 = t1072 * t1085;
t1041 = t1073 * t1084;
t1125 = t1073 * t1085;
t1123 = t1074 * t1085;
t1042 = t1075 * t1084;
t1122 = t1075 * t1085;
t1120 = t1076 * t1085;
t1043 = t1077 * t1084;
t1119 = t1077 * t1085;
t1035 = t1153 * pkin(5) - mrSges(1,2) + mrSges(3,2) + mrSges(2,3);
t1034 = t1035 * g(3);
t1037 = t1153 * pkin(1) + mrSges(1,1);
t1036 = t1037 * g(3);
t1004 = -t1034 * t1079 + (t1155 * g(3) + t1036) * t1073 + ((-t1037 - t1155) * t1079 - t1035 * t1073) * t1031;
t1118 = t1079 * t1004;
t1005 = -t1034 * t1081 + (t1156 * g(3) + t1036) * t1075 + ((-t1037 - t1156) * t1081 - t1035 * t1075) * t1032;
t1117 = t1081 * t1005;
t1006 = -t1034 * t1083 + (t1157 * g(3) + t1036) * t1077 + ((-t1037 - t1157) * t1083 - t1035 * t1077) * t1033;
t1116 = t1083 * t1006;
t1038 = pkin(1) * t1072 + qJ(3,3);
t1115 = t1085 * t1038;
t1039 = pkin(1) * t1074 + qJ(3,2);
t1114 = t1085 * t1039;
t1040 = pkin(1) * t1076 + qJ(3,1);
t1113 = t1085 * t1040;
t1109 = qJ(3,1) * t1154;
t1108 = qJ(3,2) * t1154;
t1107 = qJ(3,3) * t1154;
t1103 = (-t1078 * t1028 + t1072 * t1091) * t1136;
t1102 = (-t1080 * t1029 + t1074 * t1090) * t1135;
t1101 = (-t1082 * t1030 + t1076 * t1089) * t1134;
t1100 = t1073 * t1130;
t1099 = t1075 * t1129;
t1098 = t1077 * t1128;
t1068 = t1082 ^ 2;
t1067 = t1080 ^ 2;
t1066 = t1078 ^ 2;
t1024 = pkin(1) * qJ(3,1) - t1076 * t1128;
t1023 = pkin(1) * qJ(3,2) - t1074 * t1129;
t1022 = pkin(1) * qJ(3,3) - t1072 * t1130;
t1021 = t1095 + t1104;
t1019 = t1096 + t1105;
t1017 = t1097 + t1106;
t1015 = t1095 * t1076 + t1143;
t1014 = t1096 * t1074 + t1144;
t1013 = t1097 * t1072 + t1145;
t1 = [-g(1) * m(4) + (t1056 * t1116 + ((t1056 * t1119 - t1149) * t1068 + (t1021 * t1056 + t1053 * t1120) * t1082 + t1053 * t1040) * t1140) * t1027 + (t1055 * t1117 + ((t1055 * t1122 - t1150) * t1067 + (t1019 * t1055 + t1052 * t1123) * t1080 + t1052 * t1039) * t1141) * t1026 + (t1054 * t1118 + ((t1054 * t1125 - t1151) * t1066 + (t1017 * t1054 + t1051 * t1126) * t1078 + t1051 * t1038) * t1142) * t1025 + (-((t1053 * t1109 + t1056 * t1098) * t1068 + (-t1053 * t1024 + t1056 * t1137) * t1082 + t1015 * t1146 + t1053 * t1113) * t1101 - ((t1052 * t1108 + t1055 * t1099) * t1067 + (-t1052 * t1023 + t1055 * t1138) * t1080 + t1014 * t1147 + t1052 * t1114) * t1102 - ((t1051 * t1107 + t1054 * t1100) * t1066 + (-t1051 * t1022 + t1054 * t1139) * t1078 + t1013 * t1148 + t1051 * t1115) * t1103) * m(3); -g(2) * m(4) + (-t1053 * t1116 + ((-t1053 * t1119 - t1146) * t1068 + (-t1021 * t1053 + t1056 * t1120) * t1082 + t1056 * t1040) * t1140) * t1027 + (-t1052 * t1117 + ((-t1052 * t1122 - t1147) * t1067 + (-t1019 * t1052 + t1055 * t1123) * t1080 + t1055 * t1039) * t1141) * t1026 + (-t1051 * t1118 + ((-t1051 * t1125 - t1148) * t1066 + (-t1017 * t1051 + t1054 * t1126) * t1078 + t1054 * t1038) * t1142) * t1025 + (-((-t1053 * t1098 + t1056 * t1109) * t1068 + (-t1056 * t1024 - t1053 * t1137) * t1082 - t1015 * t1149 + t1056 * t1113) * t1101 - ((-t1052 * t1099 + t1055 * t1108) * t1067 + (-t1055 * t1023 - t1052 * t1138) * t1080 - t1014 * t1150 + t1055 * t1114) * t1102 - ((-t1051 * t1100 + t1054 * t1107) * t1066 + (-t1054 * t1022 - t1051 * t1139) * t1078 - t1013 * t1151 + t1054 * t1115) * t1103) * m(3); -t1077 * t1027 * t1006 + (t1092 * t1083 + t1043) * t1082 * t1009 * t1134 - t1075 * t1026 * t1005 + (t1093 * t1081 + t1042) * t1080 * t1008 * t1135 - t1073 * t1025 * t1004 + (t1094 * t1079 + t1041) * t1078 * t1007 * t1136 - g(3) * m(4) + (-(t1083 * t1068 * t1128 + ((0.2e1 * t1050 + pkin(1)) * t1083 + t1043) * t1110 + qJ(3,1) * (t1040 * t1083 + t1076 * t1043)) * t1101 - (t1081 * t1067 * t1129 + ((0.2e1 * t1049 + pkin(1)) * t1081 + t1042) * t1111 + qJ(3,2) * (t1039 * t1081 + t1074 * t1042)) * t1102 - (t1079 * t1066 * t1130 + ((0.2e1 * t1048 + pkin(1)) * t1079 + t1041) * t1112 + qJ(3,3) * (t1038 * t1079 + t1072 * t1041)) * t1103) * m(3);];
taugX  = t1;
