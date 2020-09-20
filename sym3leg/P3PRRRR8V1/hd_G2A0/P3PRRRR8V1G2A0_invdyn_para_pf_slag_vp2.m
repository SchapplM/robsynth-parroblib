% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G2A0
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:44
% EndTime: 2020-08-06 17:02:51
% DurationCPUTime: 5.84s
% Computational Cost: add. (18192->394), mult. (56235->755), div. (7479->7), fcn. (61926->22), ass. (0->315)
t172 = cos(qJ(3,3));
t173 = cos(qJ(2,3));
t167 = sin(qJ(2,3));
t253 = t167 * t172;
t121 = pkin(2) * t253 - pkin(5) * t173;
t155 = sin(pkin(3));
t166 = sin(qJ(3,3));
t157 = cos(pkin(3));
t338 = pkin(2) * t157;
t100 = t121 * t155 + t166 * t338;
t353 = 0.1e1 / t100;
t310 = 0.1e1 / t172 * t353;
t174 = cos(qJ(3,2));
t175 = cos(qJ(2,2));
t169 = sin(qJ(2,2));
t250 = t169 * t174;
t122 = pkin(2) * t250 - pkin(5) * t175;
t168 = sin(qJ(3,2));
t101 = t122 * t155 + t168 * t338;
t352 = 0.1e1 / t101;
t309 = 0.1e1 / t174 * t352;
t176 = cos(qJ(3,1));
t177 = cos(qJ(2,1));
t171 = sin(qJ(2,1));
t247 = t171 * t176;
t123 = pkin(2) * t247 - pkin(5) * t177;
t170 = sin(qJ(3,1));
t102 = t123 * t155 + t170 * t338;
t351 = 0.1e1 / t102;
t308 = 0.1e1 / t176 * t351;
t228 = t172 * mrSges(3,1) - mrSges(3,2) * t166;
t227 = t174 * mrSges(3,1) - mrSges(3,2) * t168;
t226 = t176 * mrSges(3,1) - mrSges(3,2) * t170;
t162 = legFrame(1,2);
t143 = sin(t162);
t146 = cos(t162);
t120 = g(1) * t146 - g(2) * t143;
t156 = cos(pkin(6));
t154 = sin(pkin(6));
t334 = g(3) * t154;
t105 = t120 * t156 - t334;
t117 = g(1) * t143 + g(2) * t146;
t217 = t105 * t157 + t117 * t155;
t161 = legFrame(2,2);
t142 = sin(t161);
t145 = cos(t161);
t119 = g(1) * t145 - g(2) * t142;
t104 = t119 * t156 - t334;
t116 = g(1) * t142 + g(2) * t145;
t218 = t104 * t157 + t116 * t155;
t160 = legFrame(3,2);
t141 = sin(t160);
t144 = cos(t160);
t118 = g(1) * t144 - g(2) * t141;
t103 = t118 * t156 - t334;
t115 = g(1) * t141 + g(2) * t144;
t219 = t103 * t157 + t115 * t155;
t354 = 2 * Ifges(3,4);
t148 = t172 ^ 2;
t350 = 0.2e1 * t148;
t150 = t174 ^ 2;
t349 = 0.2e1 * t150;
t152 = t176 ^ 2;
t348 = 0.2e1 * t152;
t347 = Ifges(3,5) / 0.2e1;
t346 = -Ifges(3,6) / 0.2e1;
t178 = xDP(3);
t183 = 0.1e1 / pkin(2);
t179 = xDP(2);
t180 = xDP(1);
t216 = t141 * t179 - t144 * t180;
t259 = t157 * t173;
t262 = t157 * t167;
t337 = pkin(2) * t172;
t79 = -(-t154 * t167 + t156 * t259) * t337 - pkin(5) * (t154 * t173 + t156 * t262);
t82 = (t154 * t259 + t156 * t167) * t337 + (t154 * t262 - t156 * t173) * pkin(5);
t58 = (t178 * t79 + t216 * t82) * t183 * t310;
t345 = pkin(2) * t58;
t215 = t142 * t179 - t145 * t180;
t258 = t157 * t175;
t261 = t157 * t169;
t336 = pkin(2) * t174;
t80 = -(-t154 * t169 + t156 * t258) * t336 - pkin(5) * (t154 * t175 + t156 * t261);
t83 = (t154 * t258 + t156 * t169) * t336 + (t154 * t261 - t156 * t175) * pkin(5);
t59 = (t178 * t80 + t215 * t83) * t183 * t309;
t344 = pkin(2) * t59;
t214 = t143 * t179 - t146 * t180;
t257 = t157 * t177;
t260 = t157 * t171;
t335 = pkin(2) * t176;
t81 = -(-t154 * t171 + t156 * t257) * t335 - pkin(5) * (t154 * t177 + t156 * t260);
t84 = (t154 * t257 + t156 * t171) * t335 + (t154 * t260 - t156 * t177) * pkin(5);
t60 = (t178 * t81 + t214 * t84) * t183 * t308;
t343 = pkin(2) * t60;
t159 = mrSges(2,2) - mrSges(3,3);
t342 = t159 / 0.2e1;
t341 = pkin(2) * t148;
t340 = pkin(2) * t150;
t339 = pkin(2) * t152;
t333 = g(3) * t156;
t181 = pkin(5) ^ 2;
t182 = pkin(2) ^ 2;
t307 = t166 * t58;
t246 = pkin(2) * t307;
t192 = t155 * t172 + t166 * t262;
t254 = t166 * t173;
t85 = t192 * t154 - t156 * t254;
t88 = -t154 * t254 - t192 * t156;
t64 = (t178 * t88 + t216 * t85) * t310;
t332 = (-pkin(5) * t246 + (t148 * t182 + t181) * t64) * t64;
t305 = t168 * t59;
t245 = pkin(2) * t305;
t191 = t155 * t174 + t168 * t261;
t251 = t168 * t175;
t86 = t191 * t154 - t156 * t251;
t89 = -t154 * t251 - t191 * t156;
t65 = (t178 * t89 + t215 * t86) * t309;
t331 = (-pkin(5) * t245 + (t150 * t182 + t181) * t65) * t65;
t303 = t170 * t60;
t244 = pkin(2) * t303;
t190 = t155 * t176 + t170 * t260;
t248 = t170 * t177;
t87 = t190 * t154 - t156 * t248;
t90 = -t154 * t248 - t190 * t156;
t66 = (t178 * t90 + t214 * t87) * t308;
t330 = (-pkin(5) * t244 + (t152 * t182 + t181) * t66) * t66;
t124 = pkin(5) * t167 + t173 * t337;
t268 = t155 * t166;
t198 = pkin(2) * t268 - t121 * t157;
t73 = -t124 * t154 + t198 * t156;
t67 = t100 * t144 + t141 * t73;
t329 = t67 * t353;
t68 = t100 * t141 - t144 * t73;
t328 = t68 * t353;
t125 = pkin(5) * t169 + t175 * t336;
t267 = t155 * t168;
t197 = pkin(2) * t267 - t122 * t157;
t74 = -t125 * t154 + t197 * t156;
t69 = t101 * t145 + t142 * t74;
t327 = t69 * t352;
t70 = t101 * t142 - t145 * t74;
t326 = t70 * t352;
t126 = pkin(5) * t171 + t177 * t335;
t266 = t155 * t170;
t196 = pkin(2) * t266 - t123 * t157;
t75 = -t126 * t154 + t196 * t156;
t71 = t102 * t146 + t143 * t75;
t325 = t71 * t351;
t72 = t102 * t143 - t146 * t75;
t324 = t72 * t351;
t76 = t124 * t156 + t198 * t154;
t323 = t76 * t353;
t77 = t125 * t156 + t197 * t154;
t322 = t77 * t352;
t78 = t126 * t156 + t196 * t154;
t321 = t78 * t351;
t225 = mrSges(3,1) * t166 + mrSges(3,2) * t172;
t91 = -t167 * t225 * t155 + t228 * t157;
t320 = t91 * t353;
t224 = mrSges(3,1) * t168 + mrSges(3,2) * t174;
t92 = -t169 * t224 * t155 + t227 * t157;
t319 = t92 * t352;
t223 = mrSges(3,1) * t170 + mrSges(3,2) * t176;
t93 = -t171 * t223 * t155 + t226 * t157;
t318 = t93 * t351;
t317 = Ifges(3,1) + Ifges(2,3);
t147 = m(1) + m(2) + m(3);
t313 = t147 * t353;
t312 = t147 * t352;
t311 = t147 * t351;
t306 = t166 * t64;
t304 = t168 * t65;
t302 = t170 * t66;
t301 = t183 * t79;
t300 = t183 * t80;
t299 = t183 * t81;
t298 = t183 * t82;
t297 = t183 * t83;
t296 = t183 * t84;
t295 = t183 * t91;
t294 = t183 * t92;
t293 = t183 * t93;
t127 = mrSges(2,1) + t228;
t289 = (t127 * t173 - t159 * t167) * t155;
t128 = mrSges(2,1) + t227;
t288 = (t128 * t175 - t159 * t169) * t155;
t129 = mrSges(2,1) + t226;
t287 = (t129 * t177 - t159 * t171) * t155;
t285 = t115 * t157;
t283 = t116 * t157;
t281 = t117 * t157;
t130 = Ifges(3,5) * t166 + Ifges(3,6) * t172;
t280 = t130 * t183;
t131 = Ifges(3,5) * t168 + Ifges(3,6) * t174;
t279 = t131 * t183;
t132 = Ifges(3,5) * t170 + Ifges(3,6) * t176;
t278 = t132 * t183;
t277 = mrSges(3,2) * t334 * t155;
t164 = xDDP(2);
t276 = t141 * t164;
t275 = t142 * t164;
t274 = t143 * t164;
t165 = xDDP(1);
t273 = t144 * t165;
t272 = t145 * t165;
t271 = t146 * t165;
t270 = t154 * t159;
t269 = t155 * t156;
t265 = t155 * t173;
t264 = t155 * t175;
t263 = t155 * t177;
t256 = t157 * t183;
t255 = t166 * t172;
t252 = t168 * t174;
t249 = t170 * t176;
t243 = pkin(5) * t306;
t242 = pkin(5) * t304;
t241 = pkin(5) * t302;
t240 = t353 * t289;
t239 = t352 * t288;
t238 = t351 * t287;
t237 = t141 * t310;
t236 = t142 * t309;
t235 = t143 * t308;
t234 = t144 * t310;
t233 = t145 * t309;
t232 = t146 * t308;
t231 = t155 * t253;
t230 = t155 * t250;
t229 = t155 * t247;
t222 = t118 * t154 + t333;
t221 = t119 * t154 + t333;
t220 = t120 * t154 + t333;
t22 = t243 - t345;
t10 = (((t157 * t58 + t64 * t265) * t341 - (-pkin(5) * t64 + t246) * t231 + t157 * t22) * t64 - (-t58 * t265 + (-t148 * t157 + t166 * t231 + t157) * t64) * t345) * t310;
t158 = Ifges(3,1) - Ifges(3,2);
t109 = -t148 * t158 + t255 * t354 + t317;
t13 = (t256 * t332 + (-t58 * t121 * t268 + t157 * (t58 * t341 - t243)) * t58) * t310;
t133 = t159 * t333;
t16 = (t22 * t345 - t332) * t353;
t4 = -t16 * t289 - t109 * t10 - t130 * t13 + 0.2e1 * t58 * ((t158 * t306 + t58 * t347) * t172 + t307 * t346 + (t350 - 0.1e1) * t64 * Ifges(3,4)) + (t118 * t270 - t219 * t127 + t133) * t173 + (t222 * t127 + t219 * t159) * t167;
t186 = t219 * t167 + t222 * t173;
t61 = t64 ^ 2;
t7 = -t91 * t16 - t130 * t10 - Ifges(3,3) * t13 - t61 * (Ifges(3,4) * t350 + t158 * t255 - Ifges(3,4)) + ((t103 * t155 - t285) * mrSges(3,1) + t186 * mrSges(3,2)) * t172 + (t277 + (-t118 * t269 + t285) * mrSges(3,2) + t186 * mrSges(3,1)) * t166;
t213 = t7 * t298 + t85 * t4;
t23 = t242 - t344;
t11 = (((t157 * t59 + t65 * t264) * t340 - (-pkin(5) * t65 + t245) * t230 + t157 * t23) * t65 + (t59 * t264 + (t150 * t157 - t168 * t230 - t157) * t65) * t344) * t309;
t110 = -t150 * t158 + t252 * t354 + t317;
t14 = (t256 * t331 + (-t59 * t122 * t267 + t157 * (t59 * t340 - t242)) * t59) * t309;
t17 = (t23 * t344 - t331) * t352;
t5 = -t17 * t288 - t110 * t11 - t131 * t14 + 0.2e1 * t59 * ((t158 * t304 + t59 * t347) * t174 + t305 * t346 + (t349 - 0.1e1) * t65 * Ifges(3,4)) + (t119 * t270 - t218 * t128 + t133) * t175 + (t221 * t128 + t218 * t159) * t169;
t185 = t218 * t169 + t221 * t175;
t62 = t65 ^ 2;
t8 = -t92 * t17 - t131 * t11 - Ifges(3,3) * t14 - t62 * (Ifges(3,4) * t349 + t158 * t252 - Ifges(3,4)) + ((t104 * t155 - t283) * mrSges(3,1) + t185 * mrSges(3,2)) * t174 + (t277 + (-t119 * t269 + t283) * mrSges(3,2) + t185 * mrSges(3,1)) * t168;
t212 = t8 * t297 + t86 * t5;
t111 = -t152 * t158 + t249 * t354 + t317;
t24 = -t241 + t343;
t12 = (((t157 * t60 + t66 * t263) * t339 - (-pkin(5) * t66 + t244) * t229 - t157 * t24) * t66 - (-t60 * t263 + (-t152 * t157 + t170 * t229 + t157) * t66) * t343) * t308;
t15 = (t256 * t330 + (-t60 * t123 * t266 + t157 * (t60 * t339 - t241)) * t60) * t308;
t18 = (-t24 * t343 - t330) * t351;
t6 = -t18 * t287 - t111 * t12 - t132 * t15 + 0.2e1 * t60 * ((t158 * t302 + t60 * t347) * t176 + t303 * t346 + (t348 - 0.1e1) * t66 * Ifges(3,4)) + (t120 * t270 - t217 * t129 + t133) * t177 + (t220 * t129 + t217 * t159) * t171;
t184 = t217 * t171 + t220 * t177;
t63 = t66 ^ 2;
t9 = -t93 * t18 - t132 * t12 - Ifges(3,3) * t15 - t63 * (Ifges(3,4) * t348 + t158 * t249 - Ifges(3,4)) + ((t105 * t155 - t281) * mrSges(3,1) + t184 * mrSges(3,2)) * t176 + (t277 + (-t120 * t269 + t281) * mrSges(3,2) + t184 * mrSges(3,1)) * t170;
t211 = t9 * t296 + t87 * t6;
t195 = t109 * t85 + t82 * t280;
t34 = t195 * t237 + t67 * t240;
t201 = Ifges(3,3) * t298 + t130 * t85;
t40 = t201 * t237 + t67 * t320;
t210 = t40 * t298 + t34 * t85;
t35 = -t195 * t234 + t68 * t240;
t41 = -t201 * t234 + t68 * t320;
t209 = t41 * t298 + t35 * t85;
t194 = t110 * t86 + t83 * t279;
t36 = t194 * t236 + t69 * t239;
t200 = Ifges(3,3) * t297 + t131 * t86;
t42 = t200 * t236 + t69 * t319;
t208 = t42 * t297 + t36 * t86;
t37 = -t194 * t233 + t70 * t239;
t43 = -t200 * t233 + t70 * t319;
t207 = t43 * t297 + t37 * t86;
t193 = t111 * t87 + t84 * t278;
t38 = t193 * t235 + t71 * t238;
t199 = Ifges(3,3) * t296 + t132 * t87;
t44 = t199 * t235 + t71 * t318;
t206 = t44 * t296 + t38 * t87;
t39 = -t193 * t232 + t72 * t238;
t45 = -t199 * t232 + t72 * t318;
t205 = t45 * t296 + t39 * t87;
t189 = t85 * t289 + t82 * t295;
t188 = t86 * t288 + t83 * t294;
t187 = t87 * t287 + t84 * t293;
t163 = xDDP(3);
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = t78 * t318 + (Ifges(3,3) * t299 + t132 * t90) * t308;
t53 = t77 * t319 + (Ifges(3,3) * t300 + t131 * t89) * t309;
t52 = t76 * t320 + (Ifges(3,3) * t301 + t130 * t88) * t310;
t51 = t78 * t238 + (t111 * t90 + t81 * t278) * t308;
t50 = t77 * t239 + (t110 * t89 + t80 * t279) * t309;
t49 = t76 * t240 + (t109 * t88 + t79 * t280) * t310;
t48 = t78 * t311 + (t90 * t287 + t81 * t293) * t308;
t47 = t77 * t312 + (t89 * t288 + t80 * t294) * t309;
t46 = t76 * t313 + (t88 * t289 + t79 * t295) * t310;
t33 = -t187 * t232 + t72 * t311;
t32 = t187 * t235 + t71 * t311;
t31 = -t188 * t233 + t70 * t312;
t30 = t188 * t236 + t69 * t312;
t29 = -t189 * t234 + t68 * t313;
t28 = t189 * t237 + t67 * t313;
t3 = -t12 * t287 - t93 * t15 + ((-mrSges(2,1) * t63 - t226 * (t63 + t57)) * t171 - 0.2e1 * (t223 * t60 + t66 * t342) * t66 * t177) * t155 - t157 * t57 * t223 + (-t18 - t117) * t147;
t2 = -t11 * t288 - t92 * t14 + ((-mrSges(2,1) * t62 - t227 * (t62 + t56)) * t169 - 0.2e1 * (t224 * t59 + t65 * t342) * t65 * t175) * t155 - t157 * t56 * t224 + (-t17 - t116) * t147;
t1 = -t10 * t289 - t91 * t13 + ((-mrSges(2,1) * t61 - t228 * (t61 + t55)) * t167 - 0.2e1 * (t225 * t58 + t64 * t342) * t64 * t173) * t155 - t157 * t55 * t225 + (-t16 - t115) * t147;
t19 = [t1 * t328 + t2 * t326 + t3 * t324 - g(1) * m(4) + (t29 * t329 + t31 * t327 + t33 * t325) * t164 + (t29 * t323 + t31 * t322 + t33 * t321) * t163 + (t29 * t328 + t31 * t326 + t33 * t324 + m(4)) * t165 + ((t45 * t299 + t39 * t90) * t163 + t205 * t274 + (-t205 * t165 - t211) * t146) * t308 + ((t43 * t300 + t37 * t89) * t163 + t207 * t275 + (-t207 * t165 - t212) * t145) * t309 + ((t41 * t301 + t35 * t88) * t163 + t209 * t276 + (-t209 * t165 - t213) * t144) * t310; t1 * t329 + t2 * t327 + t3 * t325 - g(2) * m(4) + (t28 * t328 + t30 * t326 + t32 * t324) * t165 + (t28 * t323 + t30 * t322 + t32 * t321) * t163 + (t28 * t329 + t30 * t327 + t32 * t325 + m(4)) * t164 + ((t44 * t299 + t38 * t90) * t163 - t206 * t271 + (t206 * t164 + t211) * t143) * t308 + ((t42 * t300 + t36 * t89) * t163 - t208 * t272 + (t208 * t164 + t212) * t142) * t309 + ((t40 * t301 + t34 * t88) * t163 - t210 * t273 + (t210 * t164 + t213) * t141) * t310; t1 * t323 + t2 * t322 + t3 * t321 - g(3) * m(4) + (t48 * t324 + t47 * t326 + t46 * t328) * t165 + (t48 * t325 + t47 * t327 + t46 * t329) * t164 + (t48 * t321 + t47 * t322 + t46 * t323 + m(4)) * t163 + ((t54 * t299 + t51 * t90) * t163 + t90 * t6 + t9 * t299 + (-t271 + t274) * (t54 * t296 + t51 * t87)) * t308 + ((t53 * t300 + t50 * t89) * t163 + t89 * t5 + t8 * t300 + (-t272 + t275) * (t53 * t297 + t50 * t86)) * t309 + ((t52 * t301 + t49 * t88) * t163 + t88 * t4 + t7 * t301 + (-t273 + t276) * (t52 * t298 + t49 * t85)) * t310;];
tauX  = t19;
