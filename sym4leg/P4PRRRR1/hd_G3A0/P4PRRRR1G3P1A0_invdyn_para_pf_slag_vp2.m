% Calculate vector of inverse dynamics forces for parallel robot
% P4PRRRR1G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% xDDP [4x1]
%   Generalized platform accelerations
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
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
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [4x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: xDDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3P1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:02:46
% EndTime: 2020-03-02 19:02:57
% DurationCPUTime: 10.94s
% Computational Cost: add. (4299->455), mult. (7841->858), div. (4160->18), fcn. (8222->26), ass. (0->307)
t181 = cos(qJ(3,1));
t151 = t181 ^ 2;
t327 = (0.2e1 * t151 - 0.1e1) * Ifges(3,4);
t179 = cos(qJ(3,2));
t147 = t179 ^ 2;
t326 = (0.2e1 * t147 - 0.1e1) * Ifges(3,4);
t177 = cos(qJ(3,3));
t143 = t177 ^ 2;
t325 = (0.2e1 * t143 - 0.1e1) * Ifges(3,4);
t159 = cos(qJ(3,4));
t135 = t159 ^ 2;
t324 = (0.2e1 * t135 - 0.1e1) * Ifges(3,4);
t323 = 0.2e1 * Ifges(3,4);
t136 = 0.1e1 / t159;
t144 = 0.1e1 / t177;
t148 = 0.1e1 / t179;
t152 = 0.1e1 / t181;
t199 = 0.1e1 / pkin(2);
t322 = t199 ^ 2;
t317 = mrSges(3,1) * g(3);
t162 = mrSges(2,2) - mrSges(3,3);
t316 = t162 / 0.2e1;
t163 = legFrame(4,2);
t124 = sin(t163);
t128 = cos(t163);
t103 = g(1) * t124 + g(2) * t128;
t157 = sin(qJ(3,4));
t310 = mrSges(3,1) * t157;
t112 = mrSges(3,2) * t159 + t310;
t185 = xDP(3);
t251 = t185 * t199;
t123 = mrSges(3,2) * t251;
t139 = m(1) + m(2) + m(3);
t158 = sin(qJ(2,4));
t134 = 0.1e1 / t158;
t138 = t136 / t135;
t156 = t185 ^ 2;
t261 = t156 * t199;
t184 = xDP(4);
t186 = xDP(2);
t187 = xDP(1);
t258 = t157 * t185;
t137 = 0.1e1 / t159 ^ 2;
t160 = cos(qJ(2,4));
t278 = t137 * t160;
t188 = xP(4);
t132 = sin(t188);
t133 = cos(t188);
t191 = koppelP(4,2);
t195 = koppelP(4,1);
t95 = t132 * t195 + t133 * t191;
t99 = -t132 * t191 + t133 * t195;
t30 = (-t258 * t278 + (-t128 * (-t184 * t95 + t187) + t124 * (t184 * t99 + t186)) * t136) * t199 * t134;
t29 = t30 ^ 2;
t293 = t159 * t29;
t21 = (-pkin(2) * t293 - t138 * t261) * t134;
t260 = t156 / pkin(2) ^ 2;
t239 = t138 * t260;
t241 = t136 * t251;
t25 = t137 * t260 + t29;
t292 = t160 * t30;
t309 = mrSges(3,1) * t159;
t208 = -mrSges(3,2) * t157 + mrSges(2,1) + t309;
t77 = -t158 * t162 + t208 * t160;
t259 = t157 * t159;
t277 = t137 * t199;
t9 = ((pkin(2) * t135 * t292 - t158 * t258) * t30 * t277 + (-t158 * t30 * t259 + t160 * t241) * t138 * t251) * t134;
t1 = -t77 * t9 - 0.2e1 * (t241 * t310 + t30 * t316 + t123) * t292 + (-t25 * t309 - mrSges(2,1) * t29 + (mrSges(3,2) * t25 - t112 * t239) * t157) * t158 + (-t21 - t103) * t139;
t315 = t1 * t134;
t172 = sin(qJ(2,3));
t140 = 0.1e1 / t172;
t146 = t144 / t143;
t178 = cos(qJ(2,3));
t235 = t144 * t251;
t171 = sin(qJ(3,3));
t256 = t171 * t185;
t257 = t171 * t177;
t145 = 0.1e1 / t177 ^ 2;
t270 = t145 * t199;
t192 = koppelP(3,2);
t196 = koppelP(3,1);
t100 = -t132 * t192 + t133 * t196;
t164 = legFrame(3,2);
t125 = sin(t164);
t129 = cos(t164);
t271 = t145 * t178;
t96 = t132 * t196 + t133 * t192;
t34 = (-t256 * t271 + (-t129 * (-t184 * t96 + t187) + t125 * (t100 * t184 + t186)) * t144) * t199 * t140;
t290 = t178 * t34;
t10 = ((pkin(2) * t143 * t290 - t172 * t256) * t34 * t270 + (-t172 * t34 * t257 + t178 * t235) * t146 * t251) * t140;
t104 = g(1) * t125 + g(2) * t129;
t308 = mrSges(3,1) * t171;
t116 = mrSges(3,2) * t177 + t308;
t31 = t34 ^ 2;
t291 = t177 * t31;
t22 = (-pkin(2) * t291 - t146 * t261) * t140;
t233 = t146 * t260;
t26 = t145 * t260 + t31;
t305 = mrSges(3,1) * t177;
t207 = -mrSges(3,2) * t171 + mrSges(2,1) + t305;
t78 = -t172 * t162 + t207 * t178;
t2 = -t78 * t10 - 0.2e1 * (t235 * t308 + t34 * t316 + t123) * t290 + (-t26 * t305 - mrSges(2,1) * t31 + (mrSges(3,2) * t26 - t116 * t233) * t171) * t172 + (-t22 - t104) * t139;
t314 = t140 * t2;
t174 = sin(qJ(2,2));
t141 = 0.1e1 / t174;
t165 = legFrame(2,2);
t126 = sin(t165);
t130 = cos(t165);
t105 = g(1) * t126 + g(2) * t130;
t150 = t148 / t147;
t180 = cos(qJ(2,2));
t232 = t148 * t251;
t173 = sin(qJ(3,2));
t254 = t173 * t185;
t255 = t173 * t179;
t149 = 0.1e1 / t179 ^ 2;
t266 = t149 * t199;
t193 = koppelP(2,2);
t197 = koppelP(2,1);
t101 = -t132 * t193 + t133 * t197;
t267 = t149 * t180;
t97 = t132 * t197 + t133 * t193;
t35 = (-t254 * t267 + (-t130 * (-t184 * t97 + t187) + t126 * (t101 * t184 + t186)) * t148) * t199 * t141;
t288 = t180 * t35;
t11 = ((pkin(2) * t147 * t288 - t174 * t254) * t35 * t266 + (-t174 * t35 * t255 + t180 * t232) * t150 * t251) * t141;
t307 = mrSges(3,1) * t173;
t117 = mrSges(3,2) * t179 + t307;
t32 = t35 ^ 2;
t289 = t179 * t32;
t23 = (-pkin(2) * t289 - t150 * t261) * t141;
t230 = t150 * t260;
t27 = t149 * t260 + t32;
t304 = mrSges(3,1) * t179;
t206 = -mrSges(3,2) * t173 + mrSges(2,1) + t304;
t79 = -t174 * t162 + t206 * t180;
t3 = -t79 * t11 - 0.2e1 * (t232 * t307 + t35 * t316 + t123) * t288 + (-t27 * t304 - mrSges(2,1) * t32 + (mrSges(3,2) * t27 - t117 * t230) * t173) * t174 + (-t23 - t105) * t139;
t313 = t141 * t3;
t176 = sin(qJ(2,1));
t142 = 0.1e1 / t176;
t166 = legFrame(1,2);
t127 = sin(t166);
t131 = cos(t166);
t106 = g(1) * t127 + g(2) * t131;
t175 = sin(qJ(3,1));
t306 = mrSges(3,1) * t175;
t118 = mrSges(3,2) * t181 + t306;
t154 = t152 / t151;
t182 = cos(qJ(2,1));
t229 = t152 * t251;
t252 = t175 * t185;
t253 = t175 * t181;
t153 = 0.1e1 / t181 ^ 2;
t262 = t153 * t199;
t194 = koppelP(1,2);
t198 = koppelP(1,1);
t102 = -t132 * t194 + t133 * t198;
t263 = t153 * t182;
t98 = t132 * t198 + t133 * t194;
t36 = (-t252 * t263 + (-t131 * (-t184 * t98 + t187) + t127 * (t102 * t184 + t186)) * t152) * t199 * t142;
t286 = t182 * t36;
t12 = ((pkin(2) * t151 * t286 - t176 * t252) * t36 * t262 + (-t176 * t36 * t253 + t182 * t229) * t154 * t251) * t142;
t227 = t154 * t260;
t33 = t36 ^ 2;
t287 = t181 * t33;
t24 = (-pkin(2) * t287 - t154 * t261) * t142;
t28 = t153 * t260 + t33;
t303 = mrSges(3,1) * t181;
t205 = -mrSges(3,2) * t175 + mrSges(2,1) + t303;
t80 = -t176 * t162 + t205 * t182;
t4 = -t80 * t12 - 0.2e1 * (t229 * t306 + t36 * t316 + t123) * t286 + (-t28 * t303 - mrSges(2,1) * t33 + (mrSges(3,2) * t28 - t118 * t227) * t175) * t176 + (-t24 - t106) * t139;
t312 = t142 * t4;
t311 = Ifges(3,1) + Ifges(2,3);
t302 = Ifges(3,3) * t199;
t155 = t184 ^ 2;
t167 = xDDP(4);
t169 = xDDP(2);
t301 = t134 * (-t155 * t95 + t167 * t99 + t169);
t170 = xDDP(1);
t300 = t134 * (-t155 * t99 - t167 * t95 + t170);
t299 = t140 * (t100 * t167 - t155 * t96 + t169);
t298 = t140 * (-t100 * t155 - t167 * t96 + t170);
t297 = t141 * (t101 * t167 - t155 * t97 + t169);
t296 = t141 * (-t101 * t155 - t167 * t97 + t170);
t295 = t142 * (t102 * t167 - t155 * t98 + t169);
t294 = t142 * (-t102 * t155 - t167 * t98 + t170);
t285 = t112 * t158;
t284 = t116 * t172;
t283 = t117 * t174;
t282 = t118 * t176;
t281 = t134 * t157;
t280 = t136 * t157;
t279 = t136 * t199;
t276 = t140 * t171;
t275 = t141 * t173;
t274 = t142 * t175;
t273 = t144 * t171;
t272 = t144 * t199;
t269 = t148 * t173;
t268 = t148 * t199;
t265 = t152 * t175;
t264 = t152 * t199;
t250 = t124 * t279;
t249 = t125 * t272;
t248 = t126 * t268;
t247 = t127 * t264;
t246 = t128 * t279;
t245 = t129 * t272;
t244 = t130 * t268;
t243 = t131 * t264;
t242 = t134 * t279;
t240 = t160 * t277;
t238 = t140 * t272;
t237 = t141 * t268;
t236 = t142 * t264;
t234 = t178 * t270;
t231 = t180 * t266;
t228 = t182 * t262;
t226 = t124 * t242;
t225 = t125 * t238;
t224 = t126 * t237;
t223 = t127 * t236;
t222 = t128 * t242;
t221 = t129 * t238;
t220 = t130 * t237;
t219 = t131 * t236;
t218 = Ifges(3,5) * t251 / 0.2e1;
t217 = -Ifges(3,6) * t251 / 0.2e1;
t216 = t240 * t281;
t215 = t234 * t276;
t214 = t231 * t275;
t213 = t228 * t274;
t107 = g(1) * t128 - g(2) * t124;
t212 = t103 * t158 + t107 * t160;
t108 = g(1) * t129 - g(2) * t125;
t211 = t104 * t172 + t108 * t178;
t109 = g(1) * t130 - g(2) * t126;
t210 = t105 * t174 + t109 * t180;
t110 = g(1) * t131 - g(2) * t127;
t209 = t106 * t176 + t110 * t182;
t190 = mrSges(4,1);
t189 = mrSges(4,2);
t183 = mrSges(3,2) * g(3);
t168 = xDDP(3);
t161 = Ifges(3,1) - Ifges(3,2);
t115 = Ifges(3,5) * t175 + Ifges(3,6) * t181;
t114 = Ifges(3,5) * t173 + Ifges(3,6) * t179;
t113 = Ifges(3,5) * t171 + Ifges(3,6) * t177;
t111 = Ifges(3,5) * t157 + Ifges(3,6) * t159;
t94 = -t132 * t189 + t133 * t190;
t93 = t132 * t190 + t133 * t189;
t92 = t127 * t176 + t131 * t182;
t91 = -t127 * t182 + t131 * t176;
t90 = t126 * t174 + t130 * t180;
t89 = -t126 * t180 + t130 * t174;
t88 = t125 * t172 + t129 * t178;
t87 = -t125 * t178 + t129 * t172;
t86 = -t151 * t161 + t253 * t323 + t311;
t85 = -t147 * t161 + t255 * t323 + t311;
t84 = -t143 * t161 + t257 * t323 + t311;
t83 = t124 * t158 + t128 * t160;
t82 = -t124 * t160 + t128 * t158;
t81 = -t135 * t161 + t259 * t323 + t311;
t68 = (t102 * t127 + t131 * t98) * t236;
t67 = (t101 * t126 + t130 * t97) * t237;
t66 = (t100 * t125 + t129 * t96) * t238;
t65 = (t124 * t99 + t128 * t95) * t242;
t64 = (t139 * t92 - t80 * t243) * t142;
t63 = (t139 * t91 + t80 * t247) * t142;
t62 = (t139 * t90 - t79 * t244) * t141;
t61 = (t139 * t89 + t79 * t248) * t141;
t60 = (t139 * t88 - t78 * t245) * t140;
t59 = (t139 * t87 + t78 * t249) * t140;
t58 = (t139 * t83 - t77 * t246) * t134;
t57 = (t139 * t82 + t77 * t250) * t134;
t56 = -t264 * t282 + (t139 * t152 - t80 * t228) * t274;
t55 = -t268 * t283 + (t139 * t148 - t79 * t231) * t275;
t54 = -t272 * t284 + (t139 * t144 - t78 * t234) * t276;
t53 = -t279 * t285 + (t136 * t139 - t77 * t240) * t281;
t52 = (t102 * t91 - t92 * t98) * t142;
t51 = (t101 * t89 - t90 * t97) * t141;
t50 = (t100 * t87 - t88 * t96) * t140;
t49 = (t82 * t99 - t83 * t95) * t134;
t48 = (-t86 * t243 + t80 * t92) * t142;
t47 = (t86 * t247 + t80 * t91) * t142;
t46 = (-t85 * t244 + t79 * t90) * t141;
t45 = (t85 * t248 + t79 * t89) * t141;
t44 = (-t84 * t245 + t78 * t88) * t140;
t43 = (t84 * t249 + t78 * t87) * t140;
t42 = (-t81 * t246 + t77 * t83) * t134;
t41 = (t81 * t250 + t77 * t82) * t134;
t40 = t115 * t264 + (t152 * t80 - t86 * t228) * t274;
t39 = t114 * t268 + (t148 * t79 - t85 * t231) * t275;
t38 = t113 * t272 + (t144 * t78 - t84 * t234) * t276;
t37 = t111 * t279 + (t136 * t77 - t81 * t240) * t281;
t20 = t139 * t52 + t68 * t80;
t19 = t139 * t51 + t67 * t79;
t18 = t139 * t50 + t66 * t78;
t17 = t139 * t49 + t65 * t77;
t16 = t52 * t80 + t68 * t86;
t15 = t51 * t79 + t67 * t85;
t14 = t50 * t78 + t66 * t84;
t13 = t49 * t77 + t65 * t81;
t8 = -t80 * t24 - t86 * t12 + t115 * t175 * t227 + 0.2e1 * ((t36 * t161 * t175 + t152 * t218) * t181 + t217 * t265 + t36 * t327) * t229 + (-t106 * t205 + t110 * t162) * t182 + t176 * (t106 * t162 + t110 * t205);
t7 = -t79 * t23 - t85 * t11 + t114 * t173 * t230 + 0.2e1 * ((t35 * t161 * t173 + t148 * t218) * t179 + t217 * t269 + t35 * t326) * t232 + (-t105 * t206 + t109 * t162) * t180 + t174 * (t105 * t162 + t109 * t206);
t6 = -t78 * t22 - t84 * t10 + t113 * t171 * t233 + 0.2e1 * ((t34 * t161 * t171 + t144 * t218) * t177 + t217 * t273 + t34 * t325) * t235 + (-t104 * t207 + t108 * t162) * t178 + t172 * (t104 * t162 + t108 * t207);
t5 = -t77 * t21 - t81 * t9 + t111 * t157 * t239 + 0.2e1 * ((t30 * t161 * t157 + t136 * t218) * t159 + t217 * t280 + t30 * t324) * t241 + (-t103 * t208 + t107 * t162) * t160 + t158 * (t103 * t162 + t107 * t208);
t69 = [(-t48 * t243 + t64 * t92) * t294 + (t48 * t247 + t64 * t91) * t295 + t92 * t312 - t8 * t219 + (-t46 * t244 + t62 * t90) * t296 + (t46 * t248 + t62 * t89) * t297 + t90 * t313 - t7 * t220 + (-t44 * t245 + t60 * t88) * t298 + (t44 * t249 + t60 * t87) * t299 + t88 * t314 - t6 * t221 + (-t42 * t246 + t58 * t83) * t300 + (t42 * t250 + t58 * t82) * t301 + t83 * t315 - t5 * t222 - t93 * t167 - t155 * t94 + (t170 - g(1)) * m(4) + ((-t115 * t219 - t118 * t92) * t264 + (t152 * t64 - t48 * t228) * t274 + (-t114 * t220 - t117 * t90) * t268 + (t148 * t62 - t46 * t231) * t275 + (-t113 * t221 - t116 * t88) * t272 + (t144 * t60 - t44 * t234) * t276 + (-t111 * t222 - t112 * t83) * t279 + (t136 * t58 - t42 * t240) * t281) * t168; (-t47 * t243 + t63 * t92) * t294 + (t47 * t247 + t63 * t91) * t295 + t91 * t312 + t8 * t223 + (-t45 * t244 + t61 * t90) * t296 + (t45 * t248 + t61 * t89) * t297 + t89 * t313 + t7 * t224 + (-t43 * t245 + t59 * t88) * t298 + (t43 * t249 + t59 * t87) * t299 + t87 * t314 + t6 * t225 + (-t41 * t246 + t57 * t83) * t300 + (t41 * t250 + t57 * t82) * t301 + t82 * t315 + t5 * t226 + t94 * t167 - t155 * t93 + (t169 - g(2)) * m(4) + ((t115 * t223 - t118 * t91) * t264 + (t152 * t63 - t47 * t228) * t274 + (t114 * t224 - t117 * t89) * t268 + (t148 * t61 - t45 * t231) * t275 + (t113 * t225 - t116 * t87) * t272 + (t144 * t59 - t43 * t234) * t276 + (t111 * t226 - t112 * t82) * t279 + (t136 * t57 - t41 * t240) * t281) * t168; (t24 * t282 - t115 * t12 + (mrSges(3,2) * t209 - t317) * t181 + (mrSges(3,1) * t209 + Ifges(3,3) * t227 - t161 * t287 + t183) * t175 - t33 * t327) * t264 + (t23 * t283 - t114 * t11 + (mrSges(3,2) * t210 - t317) * t179 + (mrSges(3,1) * t210 + Ifges(3,3) * t230 - t161 * t289 + t183) * t173 - t32 * t326) * t268 + (t22 * t284 - t113 * t10 + (mrSges(3,2) * t211 - t317) * t177 + (mrSges(3,1) * t211 + Ifges(3,3) * t233 - t161 * t291 + t183) * t171 - t31 * t325) * t272 + (t21 * t285 - t111 * t9 + (mrSges(3,2) * t212 - t317) * t159 + (mrSges(3,1) * t212 + Ifges(3,3) * t239 - t161 * t293 + t183) * t157 - t29 * t324) * t279 + (-t37 * t246 + t53 * t83) * t300 + (t37 * t250 + t53 * t82) * t301 + t265 * t312 + t269 * t313 + t273 * t314 + t280 * t315 + (-t40 * t243 + t56 * t92) * t294 + (t40 * t247 + t56 * t91) * t295 + (-t39 * t244 + t55 * t90) * t296 + (t39 * t248 + t55 * t89) * t297 + (-t38 * t245 + t54 * t88) * t298 + (t38 * t249 + t54 * t87) * t299 - t8 * t213 - t7 * t214 - t6 * t215 - t5 * t216 - g(3) * m(4) + (m(4) - t40 * t213 + ((-t118 * t175 + t302) * t264 + (-t115 * t322 * t263 + t56) * t274) * t152 - t39 * t214 + ((-t117 * t173 + t302) * t268 + (-t114 * t322 * t267 + t55) * t275) * t148 - t38 * t215 + ((-t116 * t171 + t302) * t272 + (-t113 * t322 * t271 + t54) * t276) * t144 - t37 * t216 + ((-t112 * t157 + t302) * t279 + (-t111 * t322 * t278 + t53) * t281) * t136) * t168; t65 * t5 + t66 * t6 + t67 * t7 + t68 * t8 + t51 * t3 + t52 * t4 + t49 * t1 + t50 * t2 - (-g(1) * t190 - g(2) * t189) * t132 + t133 * (g(1) * t189 - g(2) * t190) + Ifges(4,3) * t167 + t94 * t169 - t93 * t170 + (-t13 * t246 + t17 * t83) * t300 + (t13 * t250 + t17 * t82) * t301 + (-t16 * t243 + t20 * t92) * t294 + (t16 * t247 + t20 * t91) * t295 + (-t15 * t244 + t19 * t90) * t296 + (t15 * t248 + t19 * t89) * t297 + (-t14 * t245 + t18 * t88) * t298 + (t14 * t249 + t18 * t87) * t299 + ((t115 * t68 - t52 * t282) * t264 + (t152 * t20 - t16 * t228) * t274 + (t114 * t67 - t51 * t283) * t268 + (t148 * t19 - t15 * t231) * t275 + (t113 * t66 - t50 * t284) * t272 + (-t14 * t234 + t144 * t18) * t276 + (t111 * t65 - t49 * t285) * t279 + (-t13 * t240 + t136 * t17) * t281) * t168;];
tauX  = t69;
