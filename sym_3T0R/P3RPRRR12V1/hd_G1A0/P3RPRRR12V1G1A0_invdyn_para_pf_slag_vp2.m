% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR12V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:23
% EndTime: 2020-08-06 18:21:25
% DurationCPUTime: 2.46s
% Computational Cost: add. (5904->330), mult. (7281->553), div. (1476->11), fcn. (6312->18), ass. (0->218)
t186 = 2 * qJ(2,1);
t185 = 2 * qJ(2,2);
t184 = 2 * qJ(2,3);
t151 = (m(2) + m(3));
t247 = qJ(2,1) * t151;
t246 = qJ(2,2) * t151;
t245 = qJ(2,3) * t151;
t244 = -2 * Ifges(3,4);
t133 = sin(qJ(3,3));
t113 = 0.1e1 / t133;
t135 = sin(qJ(3,2));
t116 = 0.1e1 / t135;
t137 = sin(qJ(3,1));
t119 = 0.1e1 / t137;
t114 = 0.1e1 / t133 ^ 2;
t117 = 0.1e1 / t135 ^ 2;
t120 = 0.1e1 / t137 ^ 2;
t139 = cos(qJ(3,3));
t243 = 0.2e1 * t139;
t141 = cos(qJ(3,2));
t242 = 0.2e1 * t141;
t143 = cos(qJ(3,1));
t241 = 0.2e1 * t143;
t150 = pkin(1) + pkin(5);
t147 = xDP(3);
t199 = t113 * t147;
t183 = t139 * t199;
t148 = xDP(2);
t149 = xDP(1);
t126 = legFrame(3,3);
t101 = sin(t126);
t104 = cos(t126);
t111 = pkin(6) + t150;
t134 = sin(qJ(1,3));
t140 = cos(qJ(1,3));
t91 = pkin(3) * t133 + qJ(2,3);
t61 = t111 * t140 + t134 * t91;
t64 = t111 * t134 - t140 * t91;
t49 = t101 * t61 + t104 * t64;
t88 = 0.1e1 / t91;
t231 = t49 * t88;
t46 = -t101 * t64 + t104 * t61;
t234 = t46 * t88;
t226 = t148 * t231 + t149 * t234;
t19 = t183 + t226;
t67 = -t101 * t134 + t104 * t140;
t70 = t101 * t140 + t104 * t134;
t43 = (t148 * t70 + t149 * t67) * t88;
t219 = t111 * t43;
t237 = t43 * t88;
t10 = (-t183 + t19 - t219 + t226) * t237;
t115 = t113 * t114;
t122 = t139 ^ 2;
t129 = Ifges(3,1) - Ifges(3,2);
t223 = mrSges(3,1) * t133;
t177 = mrSges(3,2) * t139 + t223;
t228 = -mrSges(1,2) + mrSges(2,3);
t170 = t177 + t228 + t245;
t158 = 0.1e1 / pkin(3);
t188 = t147 * t158;
t189 = t147 ^ 2 / pkin(3) ^ 2;
t153 = qJ(2,3) ^ 2;
t160 = (pkin(1) ^ 2);
t227 = (mrSges(3,3) - mrSges(2,2));
t167 = m(3) * pkin(5) ^ 2 + 2 * mrSges(3,3) * pkin(5) + 2 * (m(3) * pkin(5) + t227) * pkin(1) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t31 = t129 * t122 + ((t153 + t160) * t151) + (mrSges(2,3) + t223) * t184 + (mrSges(3,2) * qJ(2,3) - Ifges(3,4) * t133) * t243 + t167;
t97 = mrSges(3,2) * t150 - Ifges(3,6);
t98 = mrSges(3,1) * t150 - Ifges(3,5);
t58 = t133 * t97 - t139 * t98;
t157 = pkin(3) ^ 2;
t171 = pkin(3) * t111 * t188;
t187 = -t111 ^ 2 - t157;
t200 = t113 * t139;
t204 = t19 * t111;
t7 = -((t139 * t219 + t199) * t133 + qJ(2,3) * t113 * t188) * t88 * t114 * t147 + (((-t171 * t200 + t204) * t133 + ((t122 * t157 - t153 + t187) * t133 + (t122 - 0.1e1) * pkin(3) * t184) * t43) * t113 + t204) * t237;
t73 = -g(1) * t101 + g(2) * t104;
t76 = g(1) * t104 + g(2) * t101;
t84 = -(m(2) * pkin(1)) - m(3) * t150 - t227;
t79 = mrSges(1,1) - t84;
t210 = t133 * mrSges(3,2);
t81 = mrSges(3,1) * t139 - t210;
t94 = mrSges(2,3) + t245;
t240 = (-t31 * t10 - t84 * t7 + (t134 * t79 - t170 * t140) * t76 + (-t170 * t134 - t79 * t140) * t73 + (-t115 * t139 * t58 + (t133 * t98 + t139 * t97) * t114) * t189 + (0.2e1 * (t177 + t94) * t19 + (t129 * t243 + (-t81 * t184 + (0.4e1 * t122 - 0.2e1) * Ifges(3,4)) * t113) * t188) * t43) * t88;
t195 = t116 * t147;
t181 = t141 * t195;
t127 = legFrame(2,3);
t102 = sin(t127);
t105 = cos(t127);
t136 = sin(qJ(1,2));
t142 = cos(qJ(1,2));
t92 = pkin(3) * t135 + qJ(2,2);
t62 = t111 * t142 + t136 * t92;
t65 = t111 * t136 - t142 * t92;
t50 = t102 * t62 + t105 * t65;
t89 = 0.1e1 / t92;
t230 = t50 * t89;
t47 = -t102 * t65 + t105 * t62;
t233 = t47 * t89;
t225 = t148 * t230 + t149 * t233;
t20 = t181 + t225;
t68 = -t102 * t136 + t105 * t142;
t71 = t102 * t142 + t105 * t136;
t44 = (t148 * t71 + t149 * t68) * t89;
t218 = t111 * t44;
t236 = t44 * t89;
t11 = (-t181 + t20 - t218 + t225) * t236;
t118 = t116 * t117;
t123 = t141 ^ 2;
t222 = mrSges(3,1) * t135;
t176 = mrSges(3,2) * t141 + t222;
t169 = t176 + t228 + t246;
t154 = qJ(2,2) ^ 2;
t32 = t129 * t123 + ((t154 + t160) * t151) + (mrSges(2,3) + t222) * t185 + (mrSges(3,2) * qJ(2,2) - Ifges(3,4) * t135) * t242 + t167;
t59 = t135 * t97 - t141 * t98;
t74 = -g(1) * t102 + g(2) * t105;
t77 = g(1) * t105 + g(2) * t102;
t196 = t116 * t141;
t203 = t20 * t111;
t8 = -((t141 * t218 + t195) * t135 + qJ(2,2) * t116 * t188) * t89 * t117 * t147 + (((-t171 * t196 + t203) * t135 + ((t123 * t157 - t154 + t187) * t135 + (t123 - 0.1e1) * pkin(3) * t185) * t44) * t116 + t203) * t236;
t209 = t135 * mrSges(3,2);
t82 = mrSges(3,1) * t141 - t209;
t95 = mrSges(2,3) + t246;
t239 = (-t32 * t11 - t84 * t8 + (t136 * t79 - t169 * t142) * t77 + (-t169 * t136 - t79 * t142) * t74 + (-t118 * t141 * t59 + (t135 * t98 + t141 * t97) * t117) * t189 + (0.2e1 * (t176 + t95) * t20 + (t129 * t242 + (-t82 * t185 + (0.4e1 * t123 - 0.2e1) * Ifges(3,4)) * t116) * t188) * t44) * t89;
t191 = t119 * t147;
t179 = t143 * t191;
t128 = legFrame(1,3);
t103 = sin(t128);
t106 = cos(t128);
t138 = sin(qJ(1,1));
t144 = cos(qJ(1,1));
t93 = pkin(3) * t137 + qJ(2,1);
t63 = t111 * t144 + t138 * t93;
t66 = t111 * t138 - t144 * t93;
t51 = t103 * t63 + t106 * t66;
t90 = 0.1e1 / t93;
t229 = t51 * t90;
t48 = -t103 * t66 + t106 * t63;
t232 = t48 * t90;
t224 = t148 * t229 + t149 * t232;
t21 = t179 + t224;
t69 = -t103 * t138 + t106 * t144;
t72 = t103 * t144 + t106 * t138;
t45 = (t148 * t72 + t149 * t69) * t90;
t217 = t111 * t45;
t235 = t45 * t90;
t12 = (-t179 + t21 - t217 + t224) * t235;
t121 = t119 * t120;
t124 = t143 ^ 2;
t221 = mrSges(3,1) * t137;
t175 = mrSges(3,2) * t143 + t221;
t168 = t175 + t228 + t247;
t155 = qJ(2,1) ^ 2;
t33 = t129 * t124 + ((t155 + t160) * t151) + (mrSges(2,3) + t221) * t186 + (mrSges(3,2) * qJ(2,1) - Ifges(3,4) * t137) * t241 + t167;
t60 = t137 * t97 - t143 * t98;
t75 = -g(1) * t103 + g(2) * t106;
t78 = g(1) * t106 + g(2) * t103;
t208 = t137 * mrSges(3,2);
t83 = mrSges(3,1) * t143 - t208;
t192 = t119 * t143;
t202 = t21 * t111;
t9 = -((t143 * t217 + t191) * t137 + qJ(2,1) * t119 * t188) * t90 * t120 * t147 + (((-t171 * t192 + t202) * t137 + ((t124 * t157 - t155 + t187) * t137 + (t124 - 0.1e1) * pkin(3) * t186) * t45) * t119 + t202) * t235;
t96 = mrSges(2,3) + t247;
t238 = (-t33 * t12 - t84 * t9 + (t138 * t79 - t168 * t144) * t78 + (-t168 * t138 - t79 * t144) * t75 + (-t121 * t143 * t60 + (t137 * t98 + t143 * t97) * t120) * t189 + (0.2e1 * (t175 + t96) * t21 + (t129 * t241 + (-t83 * t186 + (0.4e1 * t124 - 0.2e1) * Ifges(3,4)) * t119) * t188) * t45) * t90;
t220 = Ifges(3,3) * t158;
t131 = xDDP(2);
t216 = t131 * t88;
t215 = t131 * t89;
t214 = t131 * t90;
t132 = xDDP(1);
t213 = t132 * t88;
t212 = t132 * t89;
t211 = t132 * t90;
t207 = t158 * t88;
t206 = t158 * t89;
t205 = t158 * t90;
t130 = xDDP(3);
t201 = t113 * t130;
t198 = t113 * t158;
t197 = t116 * t130;
t194 = t116 * t158;
t193 = t119 * t130;
t190 = t119 * t158;
t182 = t115 * t189;
t180 = t118 * t189;
t178 = t121 * t189;
t174 = t134 * t76 - t140 * t73;
t173 = t136 * t77 - t74 * t142;
t172 = t138 * t78 - t144 * t75;
t146 = mrSges(3,1) * g(3);
t145 = mrSges(3,2) * g(3);
t57 = (t143 * t151 - t158 * t83) * t119;
t56 = (t141 * t151 - t158 * t82) * t116;
t55 = (t139 * t151 - t81 * t158) * t113;
t54 = (t143 * t84 - t158 * t60) * t119;
t53 = (t141 * t84 - t158 * t59) * t116;
t52 = (t139 * t84 - t158 * t58) * t113;
t42 = t45 ^ 2;
t41 = t44 ^ 2;
t40 = t43 ^ 2;
t30 = t120 * t189 + t42;
t29 = t117 * t189 + t41;
t28 = t114 * t189 + t40;
t27 = (t151 * t51 + t72 * t84) * t90;
t26 = (t151 * t50 + t71 * t84) * t89;
t25 = (t151 * t49 + t70 * t84) * t88;
t24 = (t151 * t48 + t69 * t84) * t90;
t23 = (t151 * t47 + t68 * t84) * t89;
t22 = (t151 * t46 + t67 * t84) * t88;
t18 = (t33 * t72 + t51 * t84) * t90;
t17 = (t32 * t71 + t50 * t84) * t89;
t16 = (t31 * t70 + t49 * t84) * t88;
t15 = (t33 * t69 + t48 * t84) * t90;
t14 = (t32 * t68 + t47 * t84) * t89;
t13 = (t31 * t67 + t46 * t84) * t88;
t6 = -t30 * t221 - t84 * t12 - t42 * t96 + (-t172 - t9) * t151 + (-mrSges(3,2) * t30 - t83 * t178) * t143;
t5 = -t29 * t222 - t84 * t11 - t41 * t95 + (-t173 - t8) * t151 + (-mrSges(3,2) * t29 - t82 * t180) * t141;
t4 = -t28 * t223 - t84 * t10 - t40 * t94 + (-t174 - t7) * t151 + (-mrSges(3,2) * t28 - t81 * t182) * t139;
t1 = [(t15 * t69 + t24 * t48) * t211 + (t15 * t72 + t24 * t51) * t214 + (t24 * t143 - (t48 * t83 + t60 * t69) * t205) * t193 + t69 * t238 + t6 * t232 + (t14 * t68 + t23 * t47) * t212 + (t14 * t71 + t23 * t50) * t215 + (t23 * t141 - (t47 * t82 + t59 * t68) * t206) * t197 + t68 * t239 + t5 * t233 + (t13 * t67 + t22 * t46) * t213 + (t13 * t70 + t22 * t49) * t216 + (t22 * t139 - (t46 * t81 + t58 * t67) * t207) * t201 + t67 * t240 + t4 * t234 + (-g(1) + t132) * m(4); (t18 * t69 + t27 * t48) * t211 + (t18 * t72 + t27 * t51) * t214 + (t27 * t143 - (t51 * t83 + t60 * t72) * t205) * t193 + t72 * t238 + t6 * t229 + (t17 * t68 + t26 * t47) * t212 + (t17 * t71 + t26 * t50) * t215 + (t26 * t141 - (t50 * t82 + t59 * t71) * t206) * t197 + t71 * t239 + t5 * t230 + (t16 * t67 + t25 * t46) * t213 + (t16 * t70 + t25 * t49) * t216 + (t25 * t139 - (t49 * t81 + t58 * t70) * t207) * t201 + t70 * t240 + t4 * t231 + (-g(2) + t131) * m(4); (t48 * t57 + t54 * t69) * t211 + (t51 * t57 + t54 * t72) * t214 + (t57 * t143 - (t143 * t83 - t220) * t190) * t193 + t6 * t192 - (-t60 * t12 - t83 * t9 - (-qJ(2,1) * t208 + t124 * t244 + Ifges(3,4)) * t42 + t137 * (t172 * mrSges(3,2) + t146) + (-Ifges(3,3) * t178 + t129 * t137 * t42 + t145 + (-qJ(2,1) * t42 - t172) * mrSges(3,1)) * t143) * t190 + (t47 * t56 + t53 * t68) * t212 + (t50 * t56 + t53 * t71) * t215 + (t56 * t141 - (t141 * t82 - t220) * t194) * t197 + t5 * t196 - (-t59 * t11 - t82 * t8 - (-qJ(2,2) * t209 + t123 * t244 + Ifges(3,4)) * t41 + t135 * (t173 * mrSges(3,2) + t146) + (-Ifges(3,3) * t180 + t129 * t135 * t41 + t145 + (-qJ(2,2) * t41 - t173) * mrSges(3,1)) * t141) * t194 + (t46 * t55 + t52 * t67) * t213 + (t49 * t55 + t52 * t70) * t216 + (t55 * t139 - (t139 * t81 - t220) * t198) * t201 + t4 * t200 - (-t58 * t10 - t81 * t7 - (-qJ(2,3) * t210 + t122 * t244 + Ifges(3,4)) * t40 + t133 * (t174 * mrSges(3,2) + t146) + (-Ifges(3,3) * t182 + t129 * t133 * t40 + t145 + (-qJ(2,3) * t40 - t174) * mrSges(3,1)) * t139) * t198 + (-g(3) + t130) * m(4);];
tauX  = t1;
