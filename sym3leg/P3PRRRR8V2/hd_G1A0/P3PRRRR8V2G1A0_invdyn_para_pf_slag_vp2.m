% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:25
% EndTime: 2020-08-06 17:35:31
% DurationCPUTime: 5.79s
% Computational Cost: add. (26184->425), mult. (50866->755), div. (2757->10), fcn. (53592->22), ass. (0->304)
t212 = cos(qJ(2,3));
t222 = pkin(7) + pkin(6);
t169 = t212 * t222;
t206 = sin(qJ(2,3));
t145 = pkin(2) * t206 - t169;
t195 = sin(pkin(4));
t211 = cos(qJ(3,3));
t197 = cos(pkin(4));
t205 = sin(qJ(3,3));
t294 = t197 * t205;
t305 = t195 * t206;
t191 = t211 ^ 2;
t350 = pkin(3) * t191;
t85 = 0.1e1 / ((pkin(3) * t294 + t145 * t195) * t211 + pkin(2) * t294 + t305 * t350);
t214 = cos(qJ(2,2));
t170 = t214 * t222;
t208 = sin(qJ(2,2));
t146 = pkin(2) * t208 - t170;
t213 = cos(qJ(3,2));
t207 = sin(qJ(3,2));
t292 = t197 * t207;
t303 = t195 * t208;
t192 = t213 ^ 2;
t349 = pkin(3) * t192;
t86 = 0.1e1 / ((pkin(3) * t292 + t146 * t195) * t213 + pkin(2) * t292 + t303 * t349);
t216 = cos(qJ(2,1));
t171 = t216 * t222;
t210 = sin(qJ(2,1));
t147 = pkin(2) * t210 - t171;
t215 = cos(qJ(3,1));
t209 = sin(qJ(3,1));
t290 = t197 * t209;
t301 = t195 * t210;
t193 = t215 ^ 2;
t348 = pkin(3) * t193;
t87 = 0.1e1 / ((pkin(3) * t290 + t147 * t195) * t215 + pkin(2) * t290 + t301 * t348);
t220 = xDP(2);
t221 = xDP(1);
t226 = 0.1e1 / pkin(3);
t166 = t211 * pkin(3) + pkin(2);
t130 = t206 * t166 - t169;
t142 = t166 * t294;
t300 = t195 * t211;
t91 = 0.1e1 / (t130 * t300 + t142);
t324 = t226 * t91;
t198 = legFrame(3,3);
t174 = sin(t198);
t177 = cos(t198);
t194 = sin(pkin(8));
t196 = cos(pkin(8));
t109 = -t194 * t174 + t177 * t196;
t112 = t196 * t174 + t177 * t194;
t283 = t206 * t222;
t312 = (t166 * t212 + t283) * t197;
t79 = -t130 * t109 - t112 * t312;
t82 = -t109 * t312 + t130 * t112;
t64 = (t220 * t79 + t221 * t82) * t324;
t371 = -0.2e1 * t64;
t167 = t213 * pkin(3) + pkin(2);
t131 = t208 * t167 - t170;
t143 = t167 * t292;
t298 = t195 * t213;
t92 = 0.1e1 / (t131 * t298 + t143);
t323 = t226 * t92;
t199 = legFrame(2,3);
t175 = sin(t199);
t178 = cos(t199);
t110 = -t194 * t175 + t178 * t196;
t113 = t196 * t175 + t178 * t194;
t281 = t208 * t222;
t311 = (t167 * t214 + t281) * t197;
t80 = -t131 * t110 - t113 * t311;
t83 = -t110 * t311 + t131 * t113;
t65 = (t220 * t80 + t221 * t83) * t323;
t370 = -0.2e1 * t65;
t168 = t215 * pkin(3) + pkin(2);
t132 = t210 * t168 - t171;
t144 = t168 * t290;
t296 = t195 * t215;
t93 = 0.1e1 / (t132 * t296 + t144);
t322 = t226 * t93;
t200 = legFrame(1,3);
t176 = sin(t200);
t179 = cos(t200);
t111 = -t194 * t176 + t179 * t196;
t114 = t196 * t176 + t179 * t194;
t279 = t210 * t222;
t310 = (t168 * t216 + t279) * t197;
t81 = -t132 * t111 - t114 * t310;
t84 = -t111 * t310 + t132 * t114;
t66 = (t220 * t81 + t221 * t84) * t322;
t369 = -0.2e1 * t66;
t345 = g(3) * t197;
t246 = t211 * mrSges(3,1) - mrSges(3,2) * t205;
t245 = t213 * mrSges(3,1) - mrSges(3,2) * t207;
t244 = t215 * mrSges(3,1) - mrSges(3,2) * t209;
t150 = pkin(2) * t216 + t279;
t302 = t195 * t209;
t232 = pkin(3) * t302 - t147 * t197;
t365 = t150 * t196 + t232 * t194;
t149 = pkin(2) * t214 + t281;
t304 = t195 * t207;
t233 = pkin(3) * t304 - t146 * t197;
t364 = t149 * t196 + t233 * t194;
t148 = pkin(2) * t212 + t283;
t306 = t195 * t205;
t234 = pkin(3) * t306 - t145 * t197;
t363 = t148 * t196 + t234 * t194;
t138 = -t176 * g(1) + t179 * g(2);
t141 = t179 * g(1) + t176 * g(2);
t236 = t138 * t196 - t141 * t194;
t346 = g(3) * t195;
t362 = t236 * t197 + t346;
t137 = -t175 * g(1) + t178 * g(2);
t140 = t178 * g(1) + t175 * g(2);
t238 = t137 * t196 - t140 * t194;
t361 = t238 * t197 + t346;
t136 = -t174 * g(1) + t177 * g(2);
t139 = t177 * g(1) + t174 * g(2);
t240 = t136 * t196 - t139 * t194;
t360 = t240 * t197 + t346;
t201 = Ifges(3,1) - Ifges(3,2);
t218 = mrSges(3,2) * pkin(2);
t356 = -2 * Ifges(3,4);
t359 = Ifges(3,4) + t193 * t356 + (-t201 * t209 + t218) * t215;
t358 = Ifges(3,4) + t192 * t356 + (-t201 * t207 + t218) * t213;
t357 = Ifges(3,4) + t191 * t356 + (-t201 * t205 + t218) * t211;
t219 = mrSges(3,1) * pkin(2);
t355 = pkin(3) * t64;
t354 = pkin(3) * t65;
t353 = pkin(3) * t66;
t180 = pkin(6) * mrSges(3,2) - Ifges(3,6);
t352 = -t180 / 0.2e1;
t181 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t351 = t181 / 0.2e1;
t189 = m(1) + m(2) + m(3);
t347 = g(3) * t189;
t227 = pkin(2) ^ 2;
t172 = t222 ^ 2 + t227;
t225 = pkin(3) ^ 2;
t269 = t205 * t355;
t277 = 0.2e1 * pkin(2) * pkin(3);
t293 = t197 * t206;
t73 = -t109 * t300 - (t109 * t293 + t212 * t112) * t205;
t76 = -t112 * t300 - (-t212 * t109 + t112 * t293) * t205;
t58 = (t220 * t76 + t221 * t73) * t85;
t344 = (-t222 * t269 + (t191 * t225 + t211 * t277 + t172) * t58) * t58;
t268 = t207 * t354;
t291 = t197 * t208;
t74 = -t110 * t298 - (t110 * t291 + t214 * t113) * t207;
t77 = -t113 * t298 - (-t214 * t110 + t113 * t291) * t207;
t59 = (t220 * t77 + t221 * t74) * t86;
t343 = (-t222 * t268 + (t192 * t225 + t213 * t277 + t172) * t59) * t59;
t267 = t209 * t353;
t289 = t197 * t210;
t75 = -t111 * t296 - (t111 * t289 + t216 * t114) * t209;
t78 = -t114 * t296 - (-t216 * t111 + t114 * t289) * t209;
t60 = (t220 * t78 + t221 * t75) * t87;
t342 = (-t222 * t267 + (t193 * t225 + t215 * t277 + t172) * t60) * t60;
t341 = t79 * t91;
t340 = t80 * t92;
t339 = t81 * t93;
t338 = t82 * t91;
t337 = t83 * t92;
t336 = t84 * t93;
t334 = mrSges(3,2) * t195;
t330 = t205 * mrSges(3,1);
t329 = t207 * mrSges(3,1);
t328 = t209 * mrSges(3,1);
t327 = t212 * t58;
t326 = t214 * t59;
t325 = t216 * t60;
t321 = t58 * t222;
t320 = t59 * t222;
t319 = t60 * t222;
t154 = t211 * mrSges(3,2) + t330;
t61 = t64 ^ 2;
t318 = t61 * t154;
t155 = t213 * mrSges(3,2) + t329;
t62 = t65 ^ 2;
t317 = t62 * t155;
t156 = t215 * mrSges(3,2) + t328;
t63 = t66 ^ 2;
t316 = t63 * t156;
t185 = m(3) * pkin(2) + mrSges(2,1);
t133 = t246 + t185;
t173 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t103 = t133 * t212 + t206 * t173;
t315 = t103 * t195;
t134 = t245 + t185;
t104 = t134 * t214 + t208 * t173;
t314 = t104 * t195;
t135 = t244 + t185;
t105 = t135 * t216 + t210 * t173;
t313 = t105 * t195;
t299 = t195 * t212;
t297 = t195 * t214;
t295 = t195 * t216;
t288 = t197 * t226;
t284 = t206 * t211;
t282 = t208 * t213;
t280 = t210 * t215;
t278 = mrSges(3,2) * t345;
t273 = -0.2e1 * t218;
t100 = -t154 * t305 + t246 * t197;
t106 = pkin(3) * t284 + t145;
t249 = t205 * t321;
t16 = t85 * t288 * t344 + (-t197 * t249 + (-t106 * t306 + t197 * (pkin(2) * t211 + t350)) * t64) / (t106 * t300 + t142) * t64;
t49 = t249 - t355;
t22 = (-t211 * t344 - (pkin(2) * t64 - t49 * t211) * t355) * t85;
t55 = t58 ^ 2;
t266 = -t100 * t16 - t189 * t22 + ((-t55 * t185 - t246 * (t55 + t61)) * t206 + (t154 * t371 + t58 * t173) * t327) * t195;
t101 = -t155 * t303 + t245 * t197;
t108 = pkin(3) * t282 + t146;
t248 = t207 * t320;
t17 = t86 * t288 * t343 + (-t197 * t248 + (-t108 * t304 + t197 * (pkin(2) * t213 + t349)) * t65) / (t108 * t298 + t143) * t65;
t50 = t248 - t354;
t23 = (-t213 * t343 - (pkin(2) * t65 - t50 * t213) * t354) * t86;
t56 = t59 ^ 2;
t265 = -t101 * t17 - t189 * t23 + ((-t56 * t185 - t245 * (t56 + t62)) * t208 + (t155 * t370 + t59 * t173) * t326) * t195;
t102 = -t156 * t301 + t244 * t197;
t107 = pkin(3) * t280 + t147;
t247 = t209 * t319;
t18 = t87 * t288 * t342 + (-t197 * t247 + (-t107 * t302 + t197 * (pkin(2) * t215 + t348)) * t66) / (t107 * t296 + t144) * t66;
t51 = t247 - t353;
t24 = (-t215 * t342 - (pkin(2) * t66 - t51 * t215) * t353) * t87;
t57 = t60 ^ 2;
t264 = -t102 * t18 - t189 * t24 + ((-t57 * t185 - t244 * (t57 + t63)) * t210 + (t156 * t369 + t60 * t173) * t325) * t195;
t255 = t100 * t324;
t115 = t194 * t293 - t196 * t212;
t118 = t194 * t212 + t196 * t293;
t258 = pkin(2) * t306;
t88 = t194 * t148 - t234 * t196;
t67 = -(t115 * t177 + t174 * t118) * t350 + (-t88 * t174 + t363 * t177) * t211 + t112 * t258;
t37 = t82 * t255 + (t189 * t67 + t73 * t315) * t85;
t254 = t101 * t323;
t116 = t194 * t291 - t196 * t214;
t119 = t194 * t214 + t196 * t291;
t257 = pkin(2) * t304;
t89 = t194 * t149 - t233 * t196;
t68 = -(t116 * t178 + t175 * t119) * t349 + (-t89 * t175 + t364 * t178) * t213 + t113 * t257;
t38 = t83 * t254 + (t189 * t68 + t74 * t314) * t86;
t253 = t102 * t322;
t117 = t194 * t289 - t196 * t216;
t120 = t194 * t216 + t196 * t289;
t256 = pkin(2) * t302;
t90 = t194 * t150 - t232 * t196;
t69 = -(t117 * t179 + t176 * t120) * t348 + (-t90 * t176 + t365 * t179) * t215 + t114 * t256;
t39 = t84 * t253 + (t189 * t69 + t75 * t313) * t87;
t263 = t39 + t38 + t37;
t70 = (-t174 * t115 + t118 * t177) * t350 + (t363 * t174 + t88 * t177) * t211 - t109 * t258;
t40 = t79 * t255 + (t189 * t70 + t76 * t315) * t85;
t71 = (-t175 * t116 + t119 * t178) * t349 + (t364 * t175 + t89 * t178) * t213 - t110 * t257;
t41 = t80 * t254 + (t189 * t71 + t77 * t314) * t86;
t72 = (-t176 * t117 + t120 * t179) * t348 + (t365 * t176 + t90 * t179) * t215 - t111 * t256;
t42 = t81 * t253 + (t189 * t72 + t78 * t313) * t87;
t262 = t42 + t41 + t40;
t261 = Ifges(3,3) * t324;
t260 = Ifges(3,3) * t323;
t259 = Ifges(3,3) * t322;
t121 = -t180 * t211 - t205 * t181;
t252 = t121 * t324;
t122 = -t180 * t213 - t207 * t181;
t251 = t122 * t323;
t123 = -t180 * t215 - t209 * t181;
t250 = t123 * t322;
t239 = t194 * t136 + t139 * t196;
t237 = t194 * t137 + t140 * t196;
t235 = t194 * t138 + t141 * t196;
t231 = Ifges(3,1) + Ifges(2,3) + (pkin(6) ^ 2 + t227) * m(3) + 0.2e1 * pkin(6) * mrSges(3,3);
t230 = t360 * t206 + t239 * t212;
t229 = t361 * t208 + t237 * t214;
t228 = t362 * t210 + t235 * t216;
t204 = xDDP(1);
t203 = xDDP(2);
t202 = xDDP(3);
t96 = -t201 * t193 + 0.2e1 * (Ifges(3,4) * t209 + t219) * t215 + t209 * t273 + t231;
t95 = -t201 * t192 + 0.2e1 * (Ifges(3,4) * t207 + t219) * t213 + t207 * t273 + t231;
t94 = -t201 * t191 + 0.2e1 * (Ifges(3,4) * t205 + t219) * t211 + t205 * t273 + t231;
t48 = t81 * t259 + (t102 * t72 + t123 * t78) * t87;
t47 = t80 * t260 + (t101 * t71 + t122 * t77) * t86;
t46 = t79 * t261 + (t100 * t70 + t121 * t76) * t85;
t45 = t84 * t259 + (t102 * t69 + t123 * t75) * t87;
t44 = t83 * t260 + (t101 * t68 + t122 * t74) * t86;
t43 = t82 * t261 + (t100 * t67 + t121 * t73) * t85;
t36 = t81 * t250 + (t72 * t313 + t78 * t96) * t87;
t35 = t80 * t251 + (t71 * t314 + t77 * t95) * t86;
t34 = t79 * t252 + (t70 * t315 + t76 * t94) * t85;
t33 = t84 * t250 + (t69 * t313 + t75 * t96) * t87;
t32 = t83 * t251 + (t68 * t314 + t74 * t95) * t86;
t31 = t82 * t252 + (t67 * t315 + t73 * t94) * t85;
t12 = (((t197 * t66 + t295 * t60) * t348 + ((-t267 + t319) * t210 + pkin(2) * t325) * t296 + t51 * t197) * t60 + (t66 * t295 + (t193 * t197 - t280 * t302 - t197) * t60) * t353) * t87;
t11 = (((t197 * t65 + t297 * t59) * t349 + ((-t268 + t320) * t208 + pkin(2) * t326) * t298 + t50 * t197) * t59 + (t65 * t297 + (t192 * t197 - t282 * t304 - t197) * t59) * t354) * t86;
t10 = (((t197 * t64 + t299 * t58) * t350 + ((-t269 + t321) * t206 + pkin(2) * t327) * t300 + t49 * t197) * t58 + (t64 * t299 + (t191 * t197 - t284 * t306 - t197) * t58) * t355) * t85;
t9 = -t102 * t24 - t123 * t12 - Ifges(3,3) * t18 + t57 * (pkin(2) * t328 + t359) + ((t195 * t236 - t345) * mrSges(3,1) + t228 * mrSges(3,2)) * t215 + (mrSges(3,1) * t228 - t236 * t334 + t278) * t209;
t8 = -t101 * t23 - t122 * t11 - Ifges(3,3) * t17 + t56 * (pkin(2) * t329 + t358) + ((t195 * t238 - t345) * mrSges(3,1) + t229 * mrSges(3,2)) * t213 + (mrSges(3,1) * t229 - t238 * t334 + t278) * t207;
t7 = -t100 * t22 - t121 * t10 - Ifges(3,3) * t16 + t55 * (pkin(2) * t330 + t357) + ((t195 * t240 - t345) * mrSges(3,1) + t230 * mrSges(3,2)) * t211 + (mrSges(3,1) * t230 - t240 * t334 + t278) * t205;
t6 = -t24 * t313 - t96 * t12 - t123 * t18 + ((t352 * t209 + t351 * t215) * t66 + (t219 * t209 + t359) * t60) * t369 + (-t135 * t362 - t173 * t235) * t216 + (t135 * t235 - t173 * t362) * t210;
t5 = -t23 * t314 - t95 * t11 - t122 * t17 + ((t352 * t207 + t351 * t213) * t65 + (t219 * t207 + t358) * t59) * t370 + (-t134 * t361 - t173 * t237) * t214 + (t134 * t237 - t173 * t361) * t208;
t4 = -t22 * t315 - t94 * t10 - t121 * t16 + ((t352 * t205 + t351 * t211) * t64 + (t219 * t205 + t357) * t58) * t371 + (-t133 * t360 - t173 * t239) * t212 + (t133 * t239 - t173 * t360) * t206;
t3 = -t12 * t313 - t197 * t316 + t264 - t347;
t2 = -t11 * t314 - t197 * t317 + t265 - t347;
t1 = -t10 * t315 - t197 * t318 + t266 - t347;
t13 = [-g(1) * m(4) + (t69 * t3 + t75 * t6) * t87 + (t68 * t2 + t74 * t5) * t86 + (t67 * t1 + t73 * t4) * t85 + (t336 * t9 + t337 * t8 + t338 * t7) * t226 + ((t33 * t78 + t39 * t72) * t87 + (t32 * t77 + t38 * t71) * t86 + (t31 * t76 + t37 * t70) * t85 + (t339 * t45 + t340 * t44 + t341 * t43) * t226) * t203 + t263 * t202 + (m(4) + (t33 * t75 + t39 * t69) * t87 + (t32 * t74 + t38 * t68) * t86 + (t31 * t73 + t37 * t67) * t85 + (t336 * t45 + t337 * t44 + t338 * t43) * t226) * t204; -g(2) * m(4) + (t72 * t3 + t78 * t6) * t87 + (t71 * t2 + t77 * t5) * t86 + (t70 * t1 + t76 * t4) * t85 + (t339 * t9 + t340 * t8 + t341 * t7) * t226 + ((t36 * t75 + t42 * t69) * t87 + (t35 * t74 + t41 * t68) * t86 + (t34 * t73 + t40 * t67) * t85 + (t336 * t48 + t337 * t47 + t338 * t46) * t226) * t204 + t262 * t202 + (m(4) + (t36 * t78 + t42 * t72) * t87 + (t35 * t77 + t41 * t71) * t86 + (t34 * t76 + t40 * t70) * t85 + (t339 * t48 + t340 * t47 + t341 * t46) * t226) * t203; t263 * t204 + t262 * t203 + (-t316 - t317 - t318) * t197 + (-t103 * t10 - t104 * t11 - t105 * t12) * t195 + t264 + t265 + t266 + (t202 - g(3)) * (0.3e1 * t189 + m(4));];
tauX  = t13;
