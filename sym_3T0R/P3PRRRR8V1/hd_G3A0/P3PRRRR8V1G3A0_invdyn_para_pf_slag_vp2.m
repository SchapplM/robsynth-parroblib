% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G3A0
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
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:04
% EndTime: 2020-08-06 17:16:10
% DurationCPUTime: 5.72s
% Computational Cost: add. (18192->391), mult. (56235->753), div. (7479->7), fcn. (61926->22), ass. (0->313)
t171 = cos(qJ(3,3));
t172 = cos(qJ(2,3));
t166 = sin(qJ(2,3));
t252 = t166 * t171;
t121 = pkin(2) * t252 - pkin(5) * t172;
t154 = sin(pkin(3));
t165 = sin(qJ(3,3));
t156 = cos(pkin(3));
t336 = pkin(2) * t156;
t100 = t121 * t154 + t165 * t336;
t351 = 0.1e1 / t100;
t308 = 0.1e1 / t171 * t351;
t173 = cos(qJ(3,2));
t174 = cos(qJ(2,2));
t168 = sin(qJ(2,2));
t249 = t168 * t173;
t122 = pkin(2) * t249 - pkin(5) * t174;
t167 = sin(qJ(3,2));
t101 = t122 * t154 + t167 * t336;
t350 = 0.1e1 / t101;
t307 = 0.1e1 / t173 * t350;
t175 = cos(qJ(3,1));
t176 = cos(qJ(2,1));
t170 = sin(qJ(2,1));
t246 = t170 * t175;
t123 = pkin(2) * t246 - pkin(5) * t176;
t169 = sin(qJ(3,1));
t102 = t123 * t154 + t169 * t336;
t349 = 0.1e1 / t102;
t306 = 0.1e1 / t175 * t349;
t155 = cos(pkin(6));
t331 = g(3) * t155;
t227 = t171 * mrSges(3,1) - mrSges(3,2) * t165;
t226 = t173 * mrSges(3,1) - mrSges(3,2) * t167;
t225 = t175 * mrSges(3,1) - mrSges(3,2) * t169;
t161 = legFrame(1,2);
t142 = sin(t161);
t145 = cos(t161);
t120 = g(1) * t145 - g(2) * t142;
t153 = sin(pkin(6));
t105 = -t120 * t153 - t331;
t117 = g(1) * t142 + g(2) * t145;
t216 = t105 * t156 + t117 * t154;
t160 = legFrame(2,2);
t141 = sin(t160);
t144 = cos(t160);
t119 = g(1) * t144 - g(2) * t141;
t104 = -t119 * t153 - t331;
t116 = g(1) * t141 + g(2) * t144;
t217 = t104 * t156 + t116 * t154;
t159 = legFrame(3,2);
t140 = sin(t159);
t143 = cos(t159);
t118 = g(1) * t143 - g(2) * t140;
t103 = -t118 * t153 - t331;
t115 = g(1) * t140 + g(2) * t143;
t218 = t103 * t156 + t115 * t154;
t352 = 2 * Ifges(3,4);
t147 = t171 ^ 2;
t348 = 0.2e1 * t147;
t149 = t173 ^ 2;
t347 = 0.2e1 * t149;
t151 = t175 ^ 2;
t346 = 0.2e1 * t151;
t345 = Ifges(3,5) / 0.2e1;
t344 = -Ifges(3,6) / 0.2e1;
t177 = xDP(3);
t182 = 0.1e1 / pkin(2);
t178 = xDP(2);
t179 = xDP(1);
t215 = t140 * t178 - t143 * t179;
t258 = t156 * t172;
t261 = t156 * t166;
t335 = pkin(2) * t171;
t79 = (-t153 * t166 + t155 * t258) * t335 + pkin(5) * (t153 * t172 + t155 * t261);
t82 = (t153 * t258 + t155 * t166) * t335 + (t153 * t261 - t155 * t172) * pkin(5);
t58 = (t177 * t82 + t215 * t79) * t182 * t308;
t343 = pkin(2) * t58;
t214 = t141 * t178 - t144 * t179;
t257 = t156 * t174;
t260 = t156 * t168;
t334 = pkin(2) * t173;
t80 = (-t153 * t168 + t155 * t257) * t334 + pkin(5) * (t153 * t174 + t155 * t260);
t83 = (t153 * t257 + t155 * t168) * t334 + (t153 * t260 - t155 * t174) * pkin(5);
t59 = (t177 * t83 + t214 * t80) * t182 * t307;
t342 = pkin(2) * t59;
t213 = t142 * t178 - t145 * t179;
t256 = t156 * t176;
t259 = t156 * t170;
t333 = pkin(2) * t175;
t81 = (-t153 * t170 + t155 * t256) * t333 + pkin(5) * (t153 * t176 + t155 * t259);
t84 = (t153 * t256 + t155 * t170) * t333 + (t153 * t259 - t155 * t176) * pkin(5);
t60 = (t177 * t84 + t213 * t81) * t182 * t306;
t341 = pkin(2) * t60;
t158 = mrSges(2,2) - mrSges(3,3);
t340 = t158 / 0.2e1;
t339 = pkin(2) * t147;
t338 = pkin(2) * t149;
t337 = pkin(2) * t151;
t332 = g(3) * t153;
t180 = pkin(5) ^ 2;
t181 = pkin(2) ^ 2;
t305 = t165 * t58;
t245 = pkin(2) * t305;
t191 = t154 * t171 + t165 * t261;
t253 = t165 * t172;
t85 = t153 * t191 - t155 * t253;
t88 = t153 * t253 + t155 * t191;
t64 = (t177 * t85 + t215 * t88) * t308;
t330 = (-pkin(5) * t245 + (t147 * t181 + t180) * t64) * t64;
t303 = t167 * t59;
t244 = pkin(2) * t303;
t190 = t154 * t173 + t167 * t260;
t250 = t167 * t174;
t86 = t153 * t190 - t155 * t250;
t89 = t153 * t250 + t155 * t190;
t65 = (t177 * t86 + t214 * t89) * t307;
t329 = (-pkin(5) * t244 + (t149 * t181 + t180) * t65) * t65;
t301 = t169 * t60;
t243 = pkin(2) * t301;
t189 = t154 * t175 + t169 * t259;
t247 = t169 * t176;
t87 = t153 * t189 - t155 * t247;
t90 = t153 * t247 + t155 * t189;
t66 = (t177 * t87 + t213 * t90) * t306;
t328 = (-pkin(5) * t243 + (t151 * t181 + t180) * t66) * t66;
t124 = pkin(5) * t166 + t172 * t335;
t267 = t154 * t165;
t197 = pkin(2) * t267 - t121 * t156;
t73 = t124 * t155 + t153 * t197;
t67 = t100 * t140 + t143 * t73;
t327 = t67 * t351;
t68 = t100 * t143 - t140 * t73;
t326 = t68 * t351;
t125 = pkin(5) * t168 + t174 * t334;
t266 = t154 * t167;
t196 = pkin(2) * t266 - t122 * t156;
t74 = t125 * t155 + t153 * t196;
t69 = t101 * t141 + t144 * t74;
t325 = t69 * t350;
t70 = t101 * t144 - t141 * t74;
t324 = t70 * t350;
t126 = pkin(5) * t170 + t176 * t333;
t265 = t154 * t169;
t195 = pkin(2) * t265 - t123 * t156;
t75 = t126 * t155 + t153 * t195;
t71 = t102 * t142 + t145 * t75;
t323 = t71 * t349;
t72 = t102 * t145 - t142 * t75;
t322 = t72 * t349;
t76 = -t124 * t153 + t155 * t197;
t321 = t76 * t351;
t77 = -t125 * t153 + t155 * t196;
t320 = t77 * t350;
t78 = -t126 * t153 + t155 * t195;
t319 = t78 * t349;
t224 = mrSges(3,1) * t165 + mrSges(3,2) * t171;
t91 = -t154 * t166 * t224 + t156 * t227;
t318 = t91 * t351;
t223 = mrSges(3,1) * t167 + mrSges(3,2) * t173;
t92 = -t154 * t168 * t223 + t156 * t226;
t317 = t92 * t350;
t222 = mrSges(3,1) * t169 + mrSges(3,2) * t175;
t93 = -t154 * t170 * t222 + t156 * t225;
t316 = t93 * t349;
t315 = Ifges(3,1) + Ifges(2,3);
t146 = m(1) + m(2) + m(3);
t311 = t146 * t351;
t310 = t146 * t350;
t309 = t146 * t349;
t304 = t165 * t64;
t302 = t167 * t65;
t300 = t169 * t66;
t299 = t182 * t79;
t298 = t182 * t80;
t297 = t182 * t81;
t296 = t182 * t82;
t295 = t182 * t83;
t294 = t182 * t84;
t293 = t182 * t91;
t292 = t182 * t92;
t291 = t182 * t93;
t127 = mrSges(2,1) + t227;
t287 = (t127 * t172 - t158 * t166) * t154;
t128 = mrSges(2,1) + t226;
t286 = (t128 * t174 - t158 * t168) * t154;
t129 = mrSges(2,1) + t225;
t285 = (t129 * t176 - t158 * t170) * t154;
t283 = t115 * t156;
t281 = t116 * t156;
t279 = t117 * t156;
t130 = Ifges(3,5) * t165 + Ifges(3,6) * t171;
t278 = t130 * t182;
t131 = Ifges(3,5) * t167 + Ifges(3,6) * t173;
t277 = t131 * t182;
t132 = Ifges(3,5) * t169 + Ifges(3,6) * t175;
t276 = t132 * t182;
t275 = mrSges(3,2) * t154 * t331;
t163 = xDDP(2);
t274 = t140 * t163;
t273 = t141 * t163;
t272 = t142 * t163;
t164 = xDDP(1);
t271 = t143 * t164;
t270 = t144 * t164;
t269 = t145 * t164;
t268 = t153 * t154;
t264 = t154 * t172;
t263 = t154 * t174;
t262 = t154 * t176;
t255 = t156 * t182;
t254 = t165 * t171;
t251 = t167 * t173;
t248 = t169 * t175;
t242 = pkin(5) * t304;
t241 = pkin(5) * t302;
t240 = pkin(5) * t300;
t239 = t351 * t287;
t238 = t350 * t286;
t237 = t349 * t285;
t236 = t140 * t308;
t235 = t141 * t307;
t234 = t142 * t306;
t233 = t143 * t308;
t232 = t144 * t307;
t231 = t145 * t306;
t230 = t154 * t252;
t229 = t154 * t249;
t228 = t154 * t246;
t221 = t118 * t155 - t332;
t220 = t119 * t155 - t332;
t219 = t120 * t155 - t332;
t22 = t242 - t343;
t10 = (((t156 * t58 + t264 * t64) * t339 - (-pkin(5) * t64 + t245) * t230 + t156 * t22) * t64 - (-t58 * t264 + (-t147 * t156 + t165 * t230 + t156) * t64) * t343) * t308;
t157 = Ifges(3,1) - Ifges(3,2);
t109 = -t147 * t157 + t254 * t352 + t315;
t13 = (t255 * t330 + (-t58 * t121 * t267 + t156 * (t58 * t339 - t242)) * t58) * t308;
t16 = (t22 * t343 - t330) * t351;
t4 = -t16 * t287 - t109 * t10 - t130 * t13 + 0.2e1 * ((t157 * t304 + t58 * t345) * t171 + t305 * t344 + (t348 - 0.1e1) * t64 * Ifges(3,4)) * t58 + (-t127 * t218 + t158 * t221) * t172 + (t127 * t221 + t158 * t218) * t166;
t185 = t218 * t166 + t221 * t172;
t61 = t64 ^ 2;
t7 = -t91 * t16 - t130 * t10 - Ifges(3,3) * t13 - (Ifges(3,4) * t348 + t157 * t254 - Ifges(3,4)) * t61 + ((t103 * t154 - t283) * mrSges(3,1) + t185 * mrSges(3,2)) * t171 + t165 * (t275 + (t118 * t268 + t283) * mrSges(3,2) + t185 * mrSges(3,1));
t212 = -t7 * t299 - t88 * t4;
t23 = t241 - t342;
t11 = (((t156 * t59 + t263 * t65) * t338 - (-pkin(5) * t65 + t244) * t229 + t156 * t23) * t65 + (t59 * t263 + (t149 * t156 - t167 * t229 - t156) * t65) * t342) * t307;
t110 = -t149 * t157 + t251 * t352 + t315;
t14 = (t255 * t329 + (-t59 * t122 * t266 + t156 * (t59 * t338 - t241)) * t59) * t307;
t17 = (t23 * t342 - t329) * t350;
t5 = -t17 * t286 - t110 * t11 - t131 * t14 + 0.2e1 * ((t157 * t302 + t59 * t345) * t173 + t303 * t344 + (t347 - 0.1e1) * t65 * Ifges(3,4)) * t59 + (-t128 * t217 + t158 * t220) * t174 + (t128 * t220 + t158 * t217) * t168;
t184 = t217 * t168 + t220 * t174;
t62 = t65 ^ 2;
t8 = -t92 * t17 - t131 * t11 - Ifges(3,3) * t14 - (Ifges(3,4) * t347 + t157 * t251 - Ifges(3,4)) * t62 + ((t104 * t154 - t281) * mrSges(3,1) + t184 * mrSges(3,2)) * t173 + t167 * (t275 + (t119 * t268 + t281) * mrSges(3,2) + t184 * mrSges(3,1));
t211 = -t8 * t298 - t89 * t5;
t111 = -t151 * t157 + t248 * t352 + t315;
t24 = t240 - t341;
t12 = (((t156 * t60 + t262 * t66) * t337 - (-pkin(5) * t66 + t243) * t228 + t156 * t24) * t66 - (-t60 * t262 + (-t151 * t156 + t169 * t228 + t156) * t66) * t341) * t306;
t15 = (t255 * t328 + (-t60 * t123 * t265 + t156 * (t60 * t337 - t240)) * t60) * t306;
t18 = (t24 * t341 - t328) * t349;
t6 = -t18 * t285 - t111 * t12 - t132 * t15 + 0.2e1 * ((t157 * t300 + t60 * t345) * t175 + t301 * t344 + (t346 - 0.1e1) * t66 * Ifges(3,4)) * t60 + (-t129 * t216 + t158 * t219) * t176 + (t129 * t219 + t158 * t216) * t170;
t183 = t216 * t170 + t219 * t176;
t63 = t66 ^ 2;
t9 = -t93 * t18 - t132 * t12 - Ifges(3,3) * t15 - (Ifges(3,4) * t346 + t157 * t248 - Ifges(3,4)) * t63 + ((t105 * t154 - t279) * mrSges(3,1) + t183 * mrSges(3,2)) * t175 + t169 * (t275 + (t120 * t268 + t279) * mrSges(3,2) + t183 * mrSges(3,1));
t210 = -t9 * t297 - t90 * t6;
t194 = t109 * t88 + t79 * t278;
t34 = -t194 * t233 + t239 * t67;
t200 = Ifges(3,3) * t299 + t130 * t88;
t40 = -t200 * t233 + t67 * t318;
t209 = t40 * t299 + t34 * t88;
t35 = t194 * t236 + t239 * t68;
t41 = t200 * t236 + t68 * t318;
t208 = t41 * t299 + t35 * t88;
t193 = t110 * t89 + t80 * t277;
t36 = -t193 * t232 + t238 * t69;
t199 = Ifges(3,3) * t298 + t131 * t89;
t42 = -t199 * t232 + t69 * t317;
t207 = t42 * t298 + t36 * t89;
t37 = t193 * t235 + t238 * t70;
t43 = t199 * t235 + t70 * t317;
t206 = t43 * t298 + t37 * t89;
t192 = t111 * t90 + t81 * t276;
t38 = -t192 * t231 + t237 * t71;
t198 = Ifges(3,3) * t297 + t132 * t90;
t44 = -t198 * t231 + t71 * t316;
t205 = t44 * t297 + t38 * t90;
t39 = t192 * t234 + t237 * t72;
t45 = t198 * t234 + t72 * t316;
t204 = t45 * t297 + t39 * t90;
t188 = t88 * t287 + t79 * t293;
t187 = t89 * t286 + t80 * t292;
t186 = t90 * t285 + t81 * t291;
t162 = xDDP(3);
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = t78 * t316 + (Ifges(3,3) * t294 + t132 * t87) * t306;
t53 = t77 * t317 + (Ifges(3,3) * t295 + t131 * t86) * t307;
t52 = t76 * t318 + (Ifges(3,3) * t296 + t130 * t85) * t308;
t51 = t78 * t237 + (t111 * t87 + t84 * t276) * t306;
t50 = t77 * t238 + (t110 * t86 + t83 * t277) * t307;
t49 = t76 * t239 + (t109 * t85 + t82 * t278) * t308;
t48 = t78 * t309 + (t87 * t285 + t84 * t291) * t306;
t47 = t77 * t310 + (t86 * t286 + t83 * t292) * t307;
t46 = t76 * t311 + (t85 * t287 + t82 * t293) * t308;
t33 = t186 * t234 + t72 * t309;
t32 = -t186 * t231 + t71 * t309;
t31 = t187 * t235 + t70 * t310;
t30 = -t187 * t232 + t69 * t310;
t29 = t188 * t236 + t68 * t311;
t28 = -t188 * t233 + t67 * t311;
t3 = -t12 * t285 - t93 * t15 + ((-mrSges(2,1) * t63 - t225 * (t63 + t57)) * t170 - 0.2e1 * t176 * t66 * (t222 * t60 + t66 * t340)) * t154 - t156 * t57 * t222 + (-t18 - t117) * t146;
t2 = -t11 * t286 - t92 * t14 + ((-mrSges(2,1) * t62 - t226 * (t62 + t56)) * t168 - 0.2e1 * t174 * t65 * (t223 * t59 + t65 * t340)) * t154 - t156 * t56 * t223 + (-t17 - t116) * t146;
t1 = -t10 * t287 - t91 * t13 + ((-mrSges(2,1) * t61 - t227 * (t61 + t55)) * t166 - 0.2e1 * t172 * t64 * (t224 * t58 + t64 * t340)) * t154 - t156 * t55 * t224 + (-t16 - t115) * t146;
t19 = [t1 * t327 + t2 * t325 + t3 * t323 - g(1) * m(4) + (t28 * t326 + t30 * t324 + t32 * t322) * t163 + (t28 * t321 + t30 * t320 + t32 * t319) * t162 + (t28 * t327 + t30 * t325 + t32 * t323 + m(4)) * t164 + ((t44 * t294 + t38 * t87) * t162 + t205 * t272 + (-t164 * t205 + t210) * t145) * t306 + ((t42 * t295 + t36 * t86) * t162 + t207 * t273 + (-t164 * t207 + t211) * t144) * t307 + ((t40 * t296 + t34 * t85) * t162 + t209 * t274 + (-t164 * t209 + t212) * t143) * t308; t1 * t326 + t2 * t324 + t3 * t322 - g(2) * m(4) + (t29 * t327 + t31 * t325 + t33 * t323) * t164 + (t29 * t321 + t31 * t320 + t33 * t319) * t162 + (t29 * t326 + t31 * t324 + t33 * t322 + m(4)) * t163 + ((t45 * t294 + t39 * t87) * t162 - t204 * t269 + (t163 * t204 - t210) * t142) * t306 + ((t43 * t295 + t37 * t86) * t162 - t206 * t270 + (t163 * t206 - t211) * t141) * t307 + ((t41 * t296 + t35 * t85) * t162 - t208 * t271 + (t163 * t208 - t212) * t140) * t308; t1 * t321 + t2 * t320 + t3 * t319 - g(3) * m(4) + (t48 * t323 + t47 * t325 + t46 * t327) * t164 + (t48 * t322 + t47 * t324 + t46 * t326) * t163 + (t48 * t319 + t47 * t320 + t46 * t321 + m(4)) * t162 + ((t54 * t294 + t51 * t87) * t162 + t87 * t6 + t9 * t294 + (-t269 + t272) * (t54 * t297 + t51 * t90)) * t306 + ((t53 * t295 + t50 * t86) * t162 + t86 * t5 + t8 * t295 + (-t270 + t273) * (t53 * t298 + t50 * t89)) * t307 + ((t52 * t296 + t49 * t85) * t162 + t85 * t4 + t7 * t296 + (-t271 + t274) * (t52 * t299 + t49 * t88)) * t308;];
tauX  = t19;
