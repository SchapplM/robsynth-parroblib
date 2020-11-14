% Calculate Gravitation load for parallel robot
% P3RPRRR9V1G2A0
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
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:45
% EndTime: 2020-08-06 18:51:45
% DurationCPUTime: 0.75s
% Computational Cost: add. (654->176), mult. (891->281), div. (42->10), fcn. (636->26), ass. (0->129)
t1165 = 2 * pkin(3);
t1085 = cos(pkin(7));
t1164 = 0.2e1 * t1085 ^ 2;
t1104 = m(2) + m(3);
t1088 = pkin(5) + qJ(2,1);
t1087 = pkin(5) + qJ(2,2);
t1086 = pkin(5) + qJ(2,3);
t1098 = cos(qJ(3,3));
t1081 = t1098 ^ 2;
t1163 = t1081 * pkin(3);
t1100 = cos(qJ(3,2));
t1082 = t1100 ^ 2;
t1162 = t1082 * pkin(3);
t1102 = cos(qJ(3,1));
t1083 = t1102 ^ 2;
t1161 = t1083 * pkin(3);
t1160 = t1098 * pkin(2);
t1159 = t1100 * pkin(2);
t1158 = t1102 * pkin(2);
t1089 = legFrame(3,2);
t1065 = sin(t1089);
t1068 = cos(t1089);
t1031 = t1065 * g(1) + t1068 * g(2);
t1078 = pkin(7) + qJ(3,3);
t1058 = cos(t1078);
t1052 = 0.1e1 / t1058;
t1055 = sin(t1078);
t1034 = t1068 * g(1) - t1065 * g(2);
t1093 = sin(qJ(1,3));
t1099 = cos(qJ(1,3));
t1119 = g(3) * t1099 + t1093 * t1034;
t1157 = ((-mrSges(3,1) * t1031 + t1119 * mrSges(3,2)) * t1058 + t1055 * (t1119 * mrSges(3,1) + mrSges(3,2) * t1031)) * t1052;
t1090 = legFrame(2,2);
t1066 = sin(t1090);
t1069 = cos(t1090);
t1032 = t1066 * g(1) + t1069 * g(2);
t1079 = pkin(7) + qJ(3,2);
t1059 = cos(t1079);
t1053 = 0.1e1 / t1059;
t1056 = sin(t1079);
t1035 = t1069 * g(1) - t1066 * g(2);
t1095 = sin(qJ(1,2));
t1101 = cos(qJ(1,2));
t1118 = g(3) * t1101 + t1095 * t1035;
t1156 = ((-mrSges(3,1) * t1032 + t1118 * mrSges(3,2)) * t1059 + t1056 * (t1118 * mrSges(3,1) + mrSges(3,2) * t1032)) * t1053;
t1091 = legFrame(1,2);
t1067 = sin(t1091);
t1070 = cos(t1091);
t1033 = t1067 * g(1) + t1070 * g(2);
t1080 = pkin(7) + qJ(3,1);
t1060 = cos(t1080);
t1054 = 0.1e1 / t1060;
t1057 = sin(t1080);
t1036 = t1070 * g(1) - t1067 * g(2);
t1097 = sin(qJ(1,1));
t1103 = cos(qJ(1,1));
t1117 = g(3) * t1103 + t1097 * t1036;
t1155 = ((-mrSges(3,1) * t1033 + t1117 * mrSges(3,2)) * t1060 + t1057 * (t1117 * mrSges(3,1) + mrSges(3,2) * t1033)) * t1054;
t1074 = pkin(6) + t1086;
t1062 = 0.1e1 / t1074;
t1154 = (-g(3) * t1093 + t1099 * t1034) * t1062;
t1075 = pkin(6) + t1087;
t1063 = 0.1e1 / t1075;
t1153 = (-g(3) * t1095 + t1101 * t1035) * t1063;
t1076 = pkin(6) + t1088;
t1064 = 0.1e1 / t1076;
t1152 = (-g(3) * t1097 + t1103 * t1036) * t1064;
t1092 = sin(qJ(3,3));
t1106 = pkin(2) / 0.2e1;
t1151 = (t1098 * pkin(3) + t1106) * t1092;
t1094 = sin(qJ(3,2));
t1150 = (t1100 * pkin(3) + t1106) * t1094;
t1096 = sin(qJ(3,1));
t1149 = (t1102 * pkin(3) + t1106) * t1096;
t1047 = pkin(1) * t1104 + mrSges(1,1);
t1046 = t1047 * g(3);
t1133 = mrSges(2,3) + mrSges(3,3) - mrSges(1,2);
t1111 = m(2) * qJ(2,3) + t1086 * m(3) + t1133;
t1084 = sin(pkin(7));
t1132 = -m(3) * pkin(2) - mrSges(2,1);
t1116 = (-mrSges(3,1) * t1098 + mrSges(3,2) * t1092 + t1132) * t1085 + (t1092 * mrSges(3,1) + mrSges(3,2) * t1098 + mrSges(2,2)) * t1084;
t1148 = t1062 * (t1093 * t1046 + (-t1093 * t1116 - t1111 * t1099) * g(3) + ((-t1047 + t1116) * t1099 - t1093 * t1111) * t1034);
t1112 = m(2) * qJ(2,2) + t1087 * m(3) + t1133;
t1115 = (-mrSges(3,1) * t1100 + mrSges(3,2) * t1094 + t1132) * t1085 + (t1094 * mrSges(3,1) + mrSges(3,2) * t1100 + mrSges(2,2)) * t1084;
t1147 = t1063 * (t1095 * t1046 + (-t1095 * t1115 - t1112 * t1101) * g(3) + ((-t1047 + t1115) * t1101 - t1095 * t1112) * t1035);
t1113 = m(2) * qJ(2,1) + t1088 * m(3) + t1133;
t1114 = (-mrSges(3,1) * t1102 + mrSges(3,2) * t1096 + t1132) * t1085 + (t1096 * mrSges(3,1) + mrSges(3,2) * t1102 + mrSges(2,2)) * t1084;
t1146 = t1064 * (t1097 * t1046 + (-t1097 * t1114 - t1113 * t1103) * g(3) + ((-t1047 + t1114) * t1103 - t1097 * t1113) * t1036);
t1145 = t1065 * t1093;
t1144 = t1066 * t1095;
t1143 = t1067 * t1097;
t1142 = t1068 * t1093;
t1141 = t1069 * t1095;
t1140 = t1070 * t1097;
t1139 = t1084 * t1092;
t1138 = t1084 * t1094;
t1137 = t1084 * t1096;
t1061 = pkin(1) * t1084;
t1136 = t1098 * (-t1092 * pkin(3) + t1061);
t1135 = t1100 * (-t1094 * pkin(3) + t1061);
t1134 = t1102 * (-t1096 * pkin(3) + t1061);
t1131 = t1052 * t1148;
t1130 = t1053 * t1147;
t1129 = t1054 * t1146;
t1128 = 0.1e1 / (t1085 * t1098 - t1139) * t1154;
t1127 = 0.1e1 / (t1085 * t1100 - t1138) * t1153;
t1126 = 0.1e1 / (t1085 * t1102 - t1137) * t1152;
t1125 = t1093 * t1139;
t1124 = t1095 * t1138;
t1123 = t1097 * t1137;
t1122 = pkin(1) * t1093 - t1099 * t1074;
t1121 = pkin(1) * t1095 - t1101 * t1075;
t1120 = pkin(1) * t1097 - t1103 * t1076;
t1110 = pkin(2) * t1125 + (t1125 * t1165 - t1122) * t1098;
t1109 = pkin(2) * t1124 + (t1124 * t1165 - t1121) * t1100;
t1108 = pkin(2) * t1123 + (t1123 * t1165 - t1120) * t1102;
t1107 = 1 / pkin(3);
t1105 = -pkin(3) / 0.2e1;
t1051 = t1085 * pkin(2) + pkin(1);
t1039 = t1161 + t1158 / 0.2e1 + t1105;
t1038 = t1162 + t1159 / 0.2e1 + t1105;
t1037 = t1163 + t1160 / 0.2e1 + t1105;
t1018 = pkin(1) * t1096 + (-pkin(3) + t1158 + 0.2e1 * t1161) * t1084;
t1017 = pkin(1) * t1094 + (-pkin(3) + t1159 + 0.2e1 * t1162) * t1084;
t1016 = pkin(1) * t1092 + (-pkin(3) + t1160 + 0.2e1 * t1163) * t1084;
t1015 = t1120 * t1137 + (t1083 - 0.1e1) * t1097 * pkin(3);
t1014 = t1121 * t1138 + (t1082 - 0.1e1) * t1095 * pkin(3);
t1013 = t1122 * t1139 + (t1081 - 0.1e1) * t1093 * pkin(3);
t1 = [(t1067 * t1057 + t1060 * t1140) * t1129 + (t1066 * t1056 + t1059 * t1141) * t1130 + (t1065 * t1055 + t1058 * t1142) * t1131 - g(1) * m(4) + (t1065 * t1157 + t1066 * t1156 + t1067 * t1155) * t1107 + (((t1039 * t1140 + t1067 * t1149) * t1164 + (t1067 * t1018 - t1108 * t1070) * t1085 - t1015 * t1070 + t1067 * t1134) * t1126 + ((t1038 * t1141 + t1066 * t1150) * t1164 + (t1066 * t1017 - t1109 * t1069) * t1085 - t1014 * t1069 + t1066 * t1135) * t1127 + ((t1037 * t1142 + t1065 * t1151) * t1164 + (t1065 * t1016 - t1110 * t1068) * t1085 - t1013 * t1068 + t1065 * t1136) * t1128) * t1104; (t1070 * t1057 - t1060 * t1143) * t1129 + (t1069 * t1056 - t1059 * t1144) * t1130 + (t1068 * t1055 - t1058 * t1145) * t1131 - g(2) * m(4) + (t1068 * t1157 + t1069 * t1156 + t1070 * t1155) * t1107 + (((-t1039 * t1143 + t1070 * t1149) * t1164 + (t1070 * t1018 + t1108 * t1067) * t1085 + t1015 * t1067 + t1070 * t1134) * t1126 + ((-t1038 * t1144 + t1069 * t1150) * t1164 + (t1069 * t1017 + t1109 * t1066) * t1085 + t1014 * t1066 + t1069 * t1135) * t1127 + ((-t1037 * t1145 + t1068 * t1151) * t1164 + (t1068 * t1016 + t1110 * t1065) * t1085 + t1013 * t1065 + t1068 * t1136) * t1128) * t1104; t1099 * t1148 + t1101 * t1147 + t1103 * t1146 - g(3) * m(4) + ((t1097 * t1076 + (pkin(3) * t1060 + t1051) * t1103) * t1152 + (t1095 * t1075 + (pkin(3) * t1059 + t1051) * t1101) * t1153 + (t1093 * t1074 + (pkin(3) * t1058 + t1051) * t1099) * t1154) * t1104;];
taugX  = t1;
