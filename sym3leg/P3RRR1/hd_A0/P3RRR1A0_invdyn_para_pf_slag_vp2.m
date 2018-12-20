% Calculate vector of inverse dynamics forces for parallel robot
% P3RRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (including platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 18:13
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX = P3RRR1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 18:13:10
% EndTime: 2018-12-20 18:13:12
% DurationCPUTime: 2.07s
% Computational Cost: add. (9579->230), mult. (10918->432), div. (2388->5), fcn. (11594->26), ass. (0->196)
t178 = 0.1e1 / pkin(1);
t143 = qJ(1,2) + qJ(2,2);
t127 = sin(t143);
t130 = cos(t143);
t154 = sin(qJ(1,2));
t160 = cos(qJ(1,2));
t235 = 0.1e1 / (t127 * t160 - t130 * t154);
t212 = t178 * t235;
t142 = qJ(1,3) + qJ(2,3);
t126 = sin(t142);
t129 = cos(t142);
t152 = sin(qJ(1,3));
t158 = cos(qJ(1,3));
t236 = 0.1e1 / (t126 * t158 - t129 * t152);
t213 = t178 * t236;
t237 = 2 * pkin(2);
t144 = qJ(1,1) + qJ(2,1);
t128 = sin(t144);
t131 = cos(t144);
t156 = sin(qJ(1,1));
t162 = cos(qJ(1,1));
t179 = (t128 * t156 + t131 * t162) * pkin(1);
t180 = (t127 * t154 + t130 * t160) * pkin(1);
t181 = (t126 * t152 + t129 * t158) * pkin(1);
t234 = 0.1e1 / (t128 * t162 - t131 * t156);
t233 = mrSges(2,1) * pkin(1);
t232 = mrSges(2,2) * pkin(1);
t151 = sin(qJ(2,3));
t157 = cos(qJ(2,3));
t231 = pkin(1) * (mrSges(2,1) * t151 + mrSges(2,2) * t157);
t153 = sin(qJ(2,2));
t159 = cos(qJ(2,2));
t230 = pkin(1) * (mrSges(2,1) * t153 + mrSges(2,2) * t159);
t155 = sin(qJ(2,1));
t161 = cos(qJ(2,1));
t229 = pkin(1) * (mrSges(2,1) * t155 + mrSges(2,2) * t161);
t123 = t151 * t232;
t177 = pkin(1) ^ 2;
t191 = m(2) * t177 + Ifges(1,3) + Ifges(2,3);
t209 = t157 * t233;
t100 = -0.2e1 * t123 + t191 + 0.2e1 * t209;
t117 = Ifges(2,3) - t123 + t209;
t166 = xP(3);
t139 = sin(t166);
t140 = cos(t166);
t169 = koppelP(3,2);
t172 = koppelP(3,1);
t105 = t139 * t172 + t140 * t169;
t108 = -t139 * t169 + t140 * t172;
t176 = 1 / pkin(2);
t210 = t176 * t178;
t197 = t236 * t210;
t145 = legFrame(3,3);
t132 = sin(t145);
t135 = cos(t145);
t85 = t126 * t135 + t129 * t132;
t61 = pkin(1) * (t132 * t158 + t135 * t152) + t85 * pkin(2);
t86 = -t126 * t132 + t129 * t135;
t64 = pkin(1) * (-t132 * t152 + t135 * t158) + t86 * pkin(2);
t28 = (-t105 * t64 + t108 * t61) * t197;
t46 = (-t105 * t86 + t108 * t85) * t213;
t228 = (t100 * t46 - t117 * t28) * t236;
t124 = t153 * t232;
t208 = t159 * t233;
t101 = -0.2e1 * t124 + t191 + 0.2e1 * t208;
t118 = Ifges(2,3) - t124 + t208;
t170 = koppelP(2,2);
t173 = koppelP(2,1);
t106 = t139 * t173 + t140 * t170;
t109 = -t139 * t170 + t140 * t173;
t196 = t235 * t210;
t146 = legFrame(2,3);
t133 = sin(t146);
t136 = cos(t146);
t87 = t127 * t136 + t130 * t133;
t62 = pkin(1) * (t133 * t160 + t136 * t154) + t87 * pkin(2);
t88 = -t127 * t133 + t130 * t136;
t65 = pkin(1) * (-t133 * t154 + t136 * t160) + t88 * pkin(2);
t29 = (-t106 * t65 + t109 * t62) * t196;
t47 = (-t106 * t88 + t109 * t87) * t212;
t227 = (t101 * t47 - t118 * t29) * t235;
t125 = t155 * t232;
t207 = t161 * t233;
t102 = -0.2e1 * t125 + t191 + 0.2e1 * t207;
t119 = Ifges(2,3) - t125 + t207;
t171 = koppelP(1,2);
t174 = koppelP(1,1);
t107 = t139 * t174 + t140 * t171;
t110 = -t139 * t171 + t140 * t174;
t195 = t234 * t210;
t147 = legFrame(1,3);
t134 = sin(t147);
t137 = cos(t147);
t89 = t128 * t137 + t131 * t134;
t63 = pkin(1) * (t134 * t162 + t137 * t156) + t89 * pkin(2);
t90 = -t128 * t134 + t131 * t137;
t66 = pkin(1) * (-t134 * t156 + t137 * t162) + t90 * pkin(2);
t30 = (-t107 * t66 + t110 * t63) * t195;
t211 = t178 * t234;
t48 = (-t107 * t90 + t110 * t89) * t211;
t226 = (t102 * t48 - t119 * t30) * t234;
t163 = xDP(3);
t164 = xDP(2);
t79 = t108 * t163 + t164;
t49 = t61 * t79 * t197;
t165 = xDP(1);
t82 = -t105 * t163 + t165;
t52 = t64 * t82 * t197;
t25 = -t52 - t49;
t37 = (t79 * t85 + t82 * t86) * t213;
t19 = t37 + t25;
t225 = t19 * t25;
t80 = t109 * t163 + t164;
t50 = t62 * t80 * t196;
t83 = -t106 * t163 + t165;
t53 = t65 * t83 * t196;
t26 = -t53 - t50;
t38 = (t80 * t87 + t83 * t88) * t212;
t20 = t38 + t26;
t224 = t20 * t26;
t81 = t110 * t163 + t164;
t51 = t63 * t81 * t195;
t84 = -t107 * t163 + t165;
t54 = t66 * t84 * t195;
t27 = -t54 - t51;
t39 = (t81 * t89 + t84 * t90) * t211;
t21 = t39 + t27;
t223 = t21 * t27;
t222 = t100 * t236;
t221 = t101 * t235;
t220 = t102 * t234;
t219 = t117 * t236;
t218 = t118 * t235;
t217 = t119 * t234;
t216 = t176 * t236;
t215 = t176 * t235;
t214 = t176 * t234;
t206 = t61 * t216;
t205 = t62 * t215;
t204 = t63 * t214;
t203 = t64 * t216;
t202 = t65 * t215;
t201 = t66 * t214;
t200 = t117 * t216;
t199 = t118 * t215;
t198 = t119 * t214;
t111 = -g(1) * t132 + g(2) * t135;
t114 = g(1) * t135 + g(2) * t132;
t194 = -t129 * (mrSges(2,1) * t111 - mrSges(2,2) * t114) + t126 * (mrSges(2,1) * t114 + mrSges(2,2) * t111);
t112 = -g(1) * t133 + g(2) * t136;
t115 = g(1) * t136 + g(2) * t133;
t193 = -t130 * (mrSges(2,1) * t112 - mrSges(2,2) * t115) + t127 * (mrSges(2,1) * t115 + mrSges(2,2) * t112);
t113 = -g(1) * t134 + g(2) * t137;
t116 = g(1) * t137 + g(2) * t134;
t192 = -t131 * (mrSges(2,1) * t113 - mrSges(2,2) * t116) + t128 * (mrSges(2,1) * t116 + mrSges(2,2) * t113);
t141 = t163 ^ 2;
t148 = xDDP(3);
t149 = xDDP(2);
t67 = -t105 * t141 + t108 * t148 + t149;
t150 = xDDP(1);
t70 = -t105 * t148 - t108 * t141 + t150;
t187 = t178 * (t61 * t67 + t64 * t70);
t68 = -t106 * t141 + t109 * t148 + t149;
t71 = -t106 * t148 - t109 * t141 + t150;
t186 = t178 * (t62 * t68 + t65 * t71);
t69 = -t107 * t141 + t110 * t148 + t149;
t72 = -t107 * t148 - t110 * t141 + t150;
t185 = t178 * (t63 * t69 + t66 * t72);
t184 = t178 * (t67 * t85 + t70 * t86);
t183 = t178 * (t68 * t87 + t71 * t88);
t182 = t178 * (t69 * t89 + t72 * t90);
t175 = pkin(2) ^ 2;
t168 = mrSges(3,1);
t167 = mrSges(3,2);
t138 = m(2) * pkin(1) + mrSges(1,1);
t104 = -t167 * t139 + t140 * t168;
t103 = t168 * t139 + t140 * t167;
t24 = -Ifges(2,3) * t30 + t119 * t48;
t23 = -Ifges(2,3) * t29 + t118 * t47;
t22 = -Ifges(2,3) * t28 + t117 * t46;
t18 = -t54 / 0.2e1 - t51 / 0.2e1 + t39;
t17 = -t53 / 0.2e1 - t50 / 0.2e1 + t38;
t16 = -t52 / 0.2e1 - t49 / 0.2e1 + t37;
t12 = (-pkin(2) * t223 + (-t21 * pkin(2) - t179 * t39) * t39) * t211;
t11 = ((-pkin(2) * t20 - t180 * t38) * t38 - pkin(2) * t224) * t212;
t10 = ((-pkin(2) * t19 - t181 * t37) * t37 - pkin(2) * t225) * t213;
t9 = ((t179 * t18 * t237 + t175 * t21 + t177 * t39) * t176 * t39 + (pkin(2) + t179) * t223) * t211;
t8 = ((t17 * t180 * t237 + t175 * t20 + t177 * t38) * t176 * t38 + (pkin(2) + t180) * t224) * t212;
t7 = ((t16 * t181 * t237 + t175 * t19 + t177 * t37) * t176 * t37 + (pkin(2) + t181) * t225) * t213;
t6 = t229 * t39 ^ 2 - Ifges(2,3) * t9 - t119 * t12 + t192;
t5 = t230 * t38 ^ 2 - Ifges(2,3) * t8 - t11 * t118 + t193;
t4 = t231 * t37 ^ 2 - Ifges(2,3) * t7 - t10 * t117 + t194;
t3 = -t102 * t12 - t119 * t9 - 0.2e1 * t27 * t18 * t229 + (mrSges(1,2) * t116 - t113 * t138) * t162 + (mrSges(1,2) * t113 + t116 * t138) * t156 + t192;
t2 = -t101 * t11 - t118 * t8 - 0.2e1 * t26 * t17 * t230 + (mrSges(1,2) * t115 - t112 * t138) * t160 + (mrSges(1,2) * t112 + t115 * t138) * t154 + t193;
t1 = -t100 * t10 - t117 * t7 - 0.2e1 * t25 * t16 * t231 + (mrSges(1,2) * t114 - t111 * t138) * t158 + (mrSges(1,2) * t111 + t114 * t138) * t152 + t194;
t13 = [-t103 * t148 - t141 * t104 + (t150 - g(1)) * m(3) + ((t90 * t3 + (-t198 * t66 + t220 * t90) * t182) * t234 + (t88 * t2 + (-t199 * t65 + t221 * t88) * t183) * t235 + (t86 * t1 + (-t200 * t64 + t222 * t86) * t184) * t236 + (-(t6 * t66 + (-Ifges(2,3) * t201 + t217 * t90) * t185) * t234 - (t5 * t65 + (-Ifges(2,3) * t202 + t218 * t88) * t186) * t235 - (t4 * t64 + (-Ifges(2,3) * t203 + t219 * t86) * t187) * t236) * t176) * t178; -t141 * t103 + t104 * t148 + (t149 - g(2)) * m(3) + ((t89 * t3 + (-t198 * t63 + t220 * t89) * t182) * t234 + (t87 * t2 + (-t199 * t62 + t221 * t87) * t183) * t235 + (t85 * t1 + (-t200 * t61 + t222 * t85) * t184) * t236 + (-(t6 * t63 + (-Ifges(2,3) * t204 + t217 * t89) * t185) * t234 - (t5 * t62 + (-Ifges(2,3) * t205 + t218 * t87) * t186) * t235 - (t4 * t61 + (-Ifges(2,3) * t206 + t219 * t85) * t187) * t236) * t176) * t178; t48 * t3 - t30 * t6 + t47 * t2 - t29 * t5 + t46 * t1 - t28 * t4 - t103 * t150 + t104 * t149 + Ifges(3,3) * t148 - (-g(1) * t168 - g(2) * t167) * t139 + t140 * (g(1) * t167 - g(2) * t168) + ((-t201 * t24 + t226 * t90) * t72 + (-t204 * t24 + t226 * t89) * t69 + (-t202 * t23 + t227 * t88) * t71 + (-t205 * t23 + t227 * t87) * t68 + (-t203 * t22 + t228 * t86) * t70 + (-t206 * t22 + t228 * t85) * t67) * t178;];
tauX  = t13;
