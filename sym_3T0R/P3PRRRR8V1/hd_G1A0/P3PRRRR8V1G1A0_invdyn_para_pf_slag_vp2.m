% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G1A0
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:51
% EndTime: 2020-08-06 16:49:56
% DurationCPUTime: 4.96s
% Computational Cost: add. (15450->367), mult. (40474->710), div. (4053->8), fcn. (48420->22), ass. (0->287)
t186 = cos(qJ(3,2));
t163 = 0.1e1 / t186;
t187 = cos(qJ(2,2));
t181 = sin(qJ(2,2));
t248 = t181 * t186;
t133 = pkin(2) * t248 - pkin(5) * t187;
t167 = sin(pkin(3));
t180 = sin(qJ(3,2));
t169 = cos(pkin(3));
t312 = pkin(2) * t169;
t212 = 0.1e1 / (t133 * t167 + t180 * t312);
t291 = t163 * t212;
t188 = cos(qJ(3,1));
t165 = 0.1e1 / t188;
t189 = cos(qJ(2,1));
t183 = sin(qJ(2,1));
t246 = t183 * t188;
t134 = pkin(2) * t246 - pkin(5) * t189;
t182 = sin(qJ(3,1));
t211 = 0.1e1 / (t134 * t167 + t182 * t312);
t279 = t211 * t165;
t306 = g(3) * t169;
t178 = sin(qJ(3,3));
t184 = cos(qJ(3,3));
t224 = t184 * mrSges(3,1) - mrSges(3,2) * t178;
t223 = t186 * mrSges(3,1) - mrSges(3,2) * t180;
t222 = t188 * mrSges(3,1) - mrSges(3,2) * t182;
t326 = 2 * Ifges(3,4);
t179 = sin(qJ(2,3));
t250 = t179 * t184;
t227 = t167 * t250;
t259 = t169 * t178;
t185 = cos(qJ(2,3));
t265 = t167 * t185;
t325 = 0.1e1 / (-pkin(5) * t265 + (t227 + t259) * pkin(2));
t160 = t184 ^ 2;
t324 = 0.2e1 * t160;
t162 = t186 ^ 2;
t323 = 0.2e1 * t162;
t164 = t188 ^ 2;
t322 = 0.2e1 * t164;
t321 = Ifges(3,5) / 0.2e1;
t320 = -Ifges(3,6) / 0.2e1;
t191 = xDP(2);
t192 = xDP(1);
t195 = 0.1e1 / pkin(2);
t161 = 0.1e1 / t184;
t294 = t161 * t325;
t233 = t195 * t294;
t170 = legFrame(3,3);
t150 = sin(t170);
t153 = cos(t170);
t166 = sin(pkin(6));
t168 = cos(pkin(6));
t108 = -t150 * t166 + t153 * t168;
t111 = t150 * t168 + t153 * t166;
t255 = t169 * t185;
t258 = t169 * t179;
t311 = pkin(2) * t184;
t73 = -(t108 * t255 - t111 * t179) * t311 - pkin(5) * (t108 * t258 + t111 * t185);
t74 = -(t108 * t179 + t111 * t255) * t311 - (-t108 * t185 + t111 * t258) * pkin(5);
t58 = (t191 * t74 + t192 * t73) * t233;
t319 = pkin(2) * t58;
t232 = t195 * t291;
t171 = legFrame(2,3);
t151 = sin(t171);
t154 = cos(t171);
t109 = -t151 * t166 + t154 * t168;
t112 = t151 * t168 + t154 * t166;
t254 = t169 * t187;
t257 = t169 * t181;
t310 = pkin(2) * t186;
t75 = -(t109 * t254 - t112 * t181) * t310 - pkin(5) * (t109 * t257 + t112 * t187);
t76 = -(t109 * t181 + t112 * t254) * t310 - (-t109 * t187 + t112 * t257) * pkin(5);
t59 = (t191 * t76 + t192 * t75) * t232;
t318 = pkin(2) * t59;
t231 = t195 * t279;
t172 = legFrame(1,3);
t152 = sin(t172);
t155 = cos(t172);
t110 = -t152 * t166 + t155 * t168;
t113 = t152 * t168 + t155 * t166;
t253 = t169 * t189;
t256 = t169 * t183;
t309 = pkin(2) * t188;
t77 = -(t110 * t253 - t113 * t183) * t309 - pkin(5) * (t110 * t256 + t113 * t189);
t78 = -(t110 * t183 + t113 * t253) * t309 - (-t110 * t189 + t113 * t256) * pkin(5);
t60 = (t191 * t78 + t192 * t77) * t231;
t317 = pkin(2) * t60;
t174 = mrSges(2,2) - mrSges(3,3);
t316 = t174 / 0.2e1;
t315 = pkin(2) * t160;
t314 = pkin(2) * t162;
t313 = pkin(2) * t164;
t159 = m(1) + m(2) + m(3);
t308 = g(3) * t159;
t307 = g(3) * t167;
t193 = pkin(5) ^ 2;
t194 = pkin(2) ^ 2;
t286 = t180 * t59;
t243 = pkin(2) * t286;
t115 = t166 * t187 + t168 * t257;
t118 = -t166 * t257 + t168 * t187;
t264 = t167 * t186;
t80 = (-t115 * t154 - t118 * t151) * t180 - t109 * t264;
t83 = (-t115 * t151 + t118 * t154) * t180 - t112 * t264;
t65 = (t191 * t83 + t192 * t80) * t291;
t305 = (-pkin(5) * t243 + (t162 * t194 + t193) * t65) * t65;
t284 = t182 * t60;
t242 = pkin(2) * t284;
t116 = t166 * t189 + t168 * t256;
t119 = -t166 * t256 + t168 * t189;
t262 = t167 * t188;
t81 = (-t116 * t155 - t119 * t152) * t182 - t110 * t262;
t84 = (-t116 * t152 + t119 * t155) * t182 - t113 * t262;
t66 = (t191 * t84 + t192 * t81) * t279;
t304 = (-pkin(5) * t242 + (t164 * t194 + t193) * t66) * t66;
t132 = pkin(2) * t250 - pkin(5) * t185;
t102 = 0.1e1 / (pkin(2) * t259 + t132 * t167);
t114 = t166 * t185 + t168 * t258;
t117 = -t166 * t258 + t168 * t185;
t266 = t167 * t184;
t79 = (-t114 * t153 - t117 * t150) * t178 - t108 * t266;
t82 = (-t114 * t150 + t117 * t153) * t178 - t111 * t266;
t64 = (t191 * t82 + t192 * t79) * t161 * t102;
t303 = t64 * t325;
t302 = Ifges(3,1) + Ifges(2,3);
t300 = mrSges(3,2) * t167;
t296 = t161 * t79;
t295 = t161 * t82;
t293 = t163 * t80;
t292 = t163 * t83;
t290 = t165 * t81;
t289 = t165 * t84;
t288 = t178 * t58;
t287 = t178 * t64;
t285 = t180 * t65;
t283 = t182 * t66;
t141 = mrSges(3,1) * t178 + mrSges(3,2) * t184;
t55 = t58 ^ 2;
t282 = t55 * t141;
t142 = mrSges(3,1) * t180 + mrSges(3,2) * t186;
t56 = t59 ^ 2;
t281 = t56 * t142;
t143 = mrSges(3,1) * t182 + mrSges(3,2) * t188;
t57 = t60 ^ 2;
t280 = t57 * t143;
t138 = mrSges(2,1) + t224;
t105 = t138 * t185 - t174 * t179;
t278 = t105 * t167;
t139 = mrSges(2,1) + t223;
t106 = t139 * t187 - t174 * t181;
t277 = t106 * t167;
t140 = mrSges(2,1) + t222;
t107 = t140 * t189 - t174 * t183;
t276 = t107 * t167;
t173 = Ifges(3,1) - Ifges(3,2);
t251 = t178 * t184;
t120 = -t160 * t173 + t251 * t326 + t302;
t275 = t120 * t161;
t249 = t180 * t186;
t121 = -t162 * t173 + t249 * t326 + t302;
t274 = t121 * t163;
t247 = t182 * t188;
t122 = -t164 * t173 + t247 * t326 + t302;
t273 = t122 * t165;
t144 = Ifges(3,5) * t178 + Ifges(3,6) * t184;
t272 = t144 * t161;
t145 = Ifges(3,5) * t180 + Ifges(3,6) * t186;
t271 = t145 * t163;
t146 = Ifges(3,5) * t182 + Ifges(3,6) * t188;
t270 = t146 * t165;
t269 = t167 * t178;
t268 = t167 * t180;
t267 = t167 * t182;
t263 = t167 * t187;
t261 = t167 * t189;
t260 = t169 * t174;
t252 = t169 * t195;
t245 = mrSges(3,2) * t306;
t244 = pkin(2) * t288;
t241 = pkin(5) * t287;
t240 = pkin(5) * t285;
t239 = pkin(5) * t283;
t28 = -pkin(5) * t244 + (t160 * t194 + t193) * t64;
t19 = (t28 * t252 * t303 + (-t58 * t132 * t269 + t169 * (t58 * t315 - t241)) * t102 * t58) * t161;
t31 = t241 - t319;
t25 = (t28 * t64 - t31 * t319) * t325;
t61 = t64 ^ 2;
t94 = -t141 * t167 * t179 + t224 * t169;
t238 = ((-mrSges(2,1) * t61 - t224 * (t61 + t55)) * t179 - 0.2e1 * t64 * (t141 * t58 + t64 * t316) * t185) * t167 - t94 * t19 + t159 * t25;
t20 = (t252 * t305 + (-t59 * t133 * t268 + t169 * (t59 * t314 - t240)) * t59) * t291;
t32 = t240 - t318;
t26 = (t32 * t318 - t305) * t212;
t62 = t65 ^ 2;
t95 = -t142 * t167 * t181 + t223 * t169;
t237 = ((-mrSges(2,1) * t62 - t223 * (t62 + t56)) * t181 - 0.2e1 * t65 * (t142 * t59 + t65 * t316) * t187) * t167 - t95 * t20 - t159 * t26;
t21 = (t252 * t304 + (-t60 * t134 * t267 + t169 * (t60 * t313 - t239)) * t60) * t279;
t33 = t239 - t317;
t27 = (t33 * t317 - t304) * t211;
t63 = t66 ^ 2;
t96 = -t143 * t167 * t183 + t222 * t169;
t236 = ((-mrSges(2,1) * t63 - t222 * (t63 + t57)) * t183 - 0.2e1 * t66 * (t143 * t60 + t66 * t316) * t189) * t167 - t96 * t21 - t159 * t27;
t219 = t94 * t233;
t230 = t161 * t278;
t135 = pkin(5) * t179 + t185 * t311;
t201 = pkin(2) * t269 - t132 * t169;
t85 = t135 * t168 + t201 * t166;
t88 = t135 * t166 - t201 * t168;
t67 = -t150 * t88 + t153 * t85;
t37 = t73 * t219 + (t159 * t67 + t79 * t230) * t102;
t218 = t95 * t232;
t229 = t163 * t277;
t136 = pkin(5) * t181 + t187 * t310;
t200 = pkin(2) * t268 - t133 * t169;
t86 = t136 * t168 + t200 * t166;
t89 = t136 * t166 - t200 * t168;
t68 = -t151 * t89 + t154 * t86;
t38 = t75 * t218 + (t159 * t68 + t80 * t229) * t212;
t216 = t96 * t231;
t228 = t165 * t276;
t137 = pkin(5) * t183 + t189 * t309;
t199 = pkin(2) * t267 - t134 * t169;
t87 = t137 * t168 + t199 * t166;
t90 = t137 * t166 - t199 * t168;
t69 = -t152 * t90 + t155 * t87;
t39 = t77 * t216 + (t159 * t69 + t81 * t228) * t211;
t235 = t39 + t38 + t37;
t70 = t150 * t85 + t153 * t88;
t40 = t74 * t219 + (t159 * t70 + t82 * t230) * t102;
t71 = t151 * t86 + t154 * t89;
t41 = t76 * t218 + (t159 * t71 + t83 * t229) * t212;
t72 = t152 * t87 + t155 * t90;
t42 = t78 * t216 + (t159 * t72 + t84 * t228) * t211;
t234 = t42 + t41 + t40;
t226 = t167 * t248;
t225 = t167 * t246;
t221 = Ifges(3,3) * t233;
t220 = Ifges(3,3) * t232;
t217 = Ifges(3,3) * t231;
t215 = t144 * t233;
t214 = t145 * t232;
t213 = t146 * t231;
t126 = -g(1) * t150 + g(2) * t153;
t129 = g(1) * t153 + g(2) * t150;
t91 = t126 * t168 - t129 * t166;
t210 = -t169 * t91 - t307;
t127 = -g(1) * t151 + g(2) * t154;
t130 = g(1) * t154 + g(2) * t151;
t92 = t127 * t168 - t130 * t166;
t209 = -t169 * t92 - t307;
t128 = -g(1) * t152 + g(2) * t155;
t131 = g(1) * t155 + g(2) * t152;
t93 = t128 * t168 - t131 * t166;
t208 = -t169 * t93 - t307;
t176 = xDDP(2);
t177 = xDDP(1);
t207 = t176 * t74 + t177 * t73;
t206 = t176 * t76 + t177 * t75;
t205 = t176 * t78 + t177 * t77;
t204 = t126 * t166 + t129 * t168;
t203 = t127 * t166 + t130 * t168;
t202 = t128 * t166 + t131 * t168;
t198 = -t210 * t179 + t204 * t185;
t197 = -t209 * t181 + t203 * t187;
t196 = -t208 * t183 + t202 * t189;
t175 = xDDP(3);
t147 = t174 * t307;
t48 = t78 * t213 + (t84 * t273 + t72 * t276) * t211;
t47 = t76 * t214 + (t83 * t274 + t71 * t277) * t212;
t46 = t74 * t215 + (t82 * t275 + t70 * t278) * t102;
t45 = t77 * t213 + (t81 * t273 + t69 * t276) * t211;
t44 = t75 * t214 + (t80 * t274 + t68 * t277) * t212;
t43 = t73 * t215 + (t79 * t275 + t67 * t278) * t102;
t12 = (((t169 * t60 + t66 * t261) * t313 - (-pkin(5) * t66 + t242) * t225 + t169 * t33) * t66 - (-t60 * t261 + (-t164 * t169 + t182 * t225 + t169) * t66) * t317) * t279;
t11 = (((t169 * t59 + t65 * t263) * t314 - (-pkin(5) * t65 + t243) * t226 + t169 * t32) * t65 + (t59 * t263 + (t162 * t169 - t180 * t226 - t169) * t65) * t318) * t291;
t10 = (((t169 * t58 + t64 * t265) * t315 - (-pkin(5) * t64 + t244) * t227 + t169 * t31) * t303 + (t58 * t265 + (t160 * t169 - t178 * t227 - t169) * t64) * t325 * t319) * t161;
t9 = -t96 * t27 - t146 * t12 - Ifges(3,3) * t21 - t63 * (Ifges(3,4) * t322 + t173 * t247 - Ifges(3,4)) + ((t93 * t167 - t306) * mrSges(3,1) + t196 * mrSges(3,2)) * t188 + (t196 * mrSges(3,1) - t93 * t300 + t245) * t182;
t8 = -t95 * t26 - t145 * t11 - Ifges(3,3) * t20 - t62 * (Ifges(3,4) * t323 + t173 * t249 - Ifges(3,4)) + ((t92 * t167 - t306) * mrSges(3,1) + t197 * mrSges(3,2)) * t186 + (t197 * mrSges(3,1) - t92 * t300 + t245) * t180;
t7 = t94 * t25 - t144 * t10 - Ifges(3,3) * t19 - t61 * (Ifges(3,4) * t324 + t173 * t251 - Ifges(3,4)) + ((t91 * t167 - t306) * mrSges(3,1) + t198 * mrSges(3,2)) * t184 + (t198 * mrSges(3,1) - t91 * t300 + t245) * t178;
t6 = -t27 * t276 - t122 * t12 - t146 * t21 + 0.2e1 * t60 * ((t173 * t283 + t60 * t321) * t188 + t284 * t320 + (t322 - 0.1e1) * t66 * Ifges(3,4)) + (t208 * t140 + t202 * t174) * t189 + (t202 * t140 + t93 * t260 + t147) * t183;
t5 = -t26 * t277 - t121 * t11 - t145 * t20 + 0.2e1 * t59 * ((t173 * t285 + t59 * t321) * t186 + t286 * t320 + (t323 - 0.1e1) * t65 * Ifges(3,4)) + (t209 * t139 + t203 * t174) * t187 + (t203 * t139 + t92 * t260 + t147) * t181;
t4 = t25 * t278 - t120 * t10 - t144 * t19 + 0.2e1 * t58 * ((t173 * t287 + t58 * t321) * t184 + t288 * t320 + (t324 - 0.1e1) * t64 * Ifges(3,4)) + (t210 * t138 + t204 * t174) * t185 + (t204 * t138 + t91 * t260 + t147) * t179;
t3 = -t12 * t276 - t169 * t280 + t236 - t308;
t2 = -t11 * t277 - t169 * t281 + t237 - t308;
t1 = -t10 * t278 - t169 * t282 + t238 - t308;
t13 = [(-g(1) + t177) * m(4) + t235 * t175 + ((t45 * t290 + t39 * t69) * t177 + (t45 * t289 + t39 * t72) * t176 + t69 * t3 + t6 * t290) * t211 + ((t44 * t293 + t38 * t68) * t177 + (t44 * t292 + t38 * t71) * t176 + t68 * t2 + t5 * t293) * t212 + ((t75 * t8 + t206 * (t75 * t220 + (t80 * t271 + t68 * t95) * t212)) * t291 + (t7 * t73 + t207 * (t73 * t221 + (t79 * t272 + t67 * t94) * t102)) * t294 + (t77 * t9 + t205 * (t77 * t217 + (t81 * t270 + t69 * t96) * t211)) * t279) * t195 + ((t43 * t296 + t37 * t67) * t177 + (t43 * t295 + t37 * t70) * t176 + t67 * t1 + t4 * t296) * t102; (-g(2) + t176) * m(4) + t234 * t175 + ((t48 * t290 + t42 * t69) * t177 + (t48 * t289 + t42 * t72) * t176 + t72 * t3 + t6 * t289) * t211 + ((t47 * t293 + t41 * t68) * t177 + (t47 * t292 + t41 * t71) * t176 + t71 * t2 + t5 * t292) * t212 + ((t76 * t8 + t206 * (t76 * t220 + (t83 * t271 + t71 * t95) * t212)) * t291 + (t7 * t74 + t207 * (t74 * t221 + (t82 * t272 + t70 * t94) * t102)) * t294 + (t78 * t9 + t205 * (t78 * t217 + (t84 * t270 + t72 * t96) * t211)) * t279) * t195 + ((t46 * t296 + t40 * t67) * t177 + (t46 * t295 + t40 * t70) * t176 + t70 * t1 + t4 * t295) * t102; t235 * t177 + t234 * t176 + (-t280 - t281 - t282) * t169 + (-t105 * t10 - t106 * t11 - t107 * t12) * t167 + t236 + t237 + t238 + (-g(3) + t175) * (m(4) + 0.3e1 * t159);];
tauX  = t13;
