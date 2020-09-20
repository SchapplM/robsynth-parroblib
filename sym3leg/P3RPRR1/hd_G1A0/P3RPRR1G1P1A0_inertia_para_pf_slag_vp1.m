% Calculate inertia matrix for parallel robot
% P3RPRR1G1P1A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRR1G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:11
% EndTime: 2020-03-09 21:23:11
% DurationCPUTime: 0.32s
% Computational Cost: add. (1044->90), mult. (1178->143), div. (108->4), fcn. (456->32), ass. (0->73)
t250 = pkin(7) + qJ(3,1);
t272 = pkin(1) * sin(t250) + sin(qJ(3,1)) * pkin(2);
t294 = (-t272 * rSges(3,2) + (pkin(1) * cos(t250) + pkin(2) * cos(qJ(3,1))) * rSges(3,1)) * m(3);
t249 = pkin(7) + qJ(3,2);
t273 = pkin(1) * sin(t249) + sin(qJ(3,2)) * pkin(2);
t293 = (-t273 * rSges(3,2) + (pkin(1) * cos(t249) + cos(qJ(3,2)) * pkin(2)) * rSges(3,1)) * m(3);
t248 = pkin(7) + qJ(3,3);
t274 = pkin(1) * sin(t248) + sin(qJ(3,3)) * pkin(2);
t292 = (-t274 * rSges(3,2) + (pkin(1) * cos(t248) + cos(qJ(3,3)) * pkin(2)) * rSges(3,1)) * m(3);
t275 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t231 = t275 * m(3) + Icges(3,3);
t208 = t292 + t231;
t263 = 0.1e1 / pkin(3);
t291 = t208 * t263;
t209 = t293 + t231;
t290 = t209 * t263;
t210 = t294 + t231;
t289 = t210 * t263;
t245 = legFrame(3,3) + qJ(1,3);
t233 = pkin(7) + t245;
t228 = qJ(3,3) + t233;
t221 = sin(t228);
t211 = -pkin(1) * sin(t245) - pkin(2) * sin(t233) - pkin(3) * t221;
t217 = 0.1e1 / t274;
t288 = t211 * t217;
t246 = legFrame(2,3) + qJ(1,2);
t234 = pkin(7) + t246;
t229 = qJ(3,2) + t234;
t222 = sin(t229);
t212 = -pkin(1) * sin(t246) - pkin(2) * sin(t234) - pkin(3) * t222;
t218 = 0.1e1 / t273;
t287 = t212 * t218;
t247 = legFrame(1,3) + qJ(1,1);
t235 = pkin(7) + t247;
t230 = qJ(3,1) + t235;
t223 = sin(t230);
t213 = -pkin(1) * sin(t247) - pkin(2) * sin(t235) - pkin(3) * t223;
t219 = 0.1e1 / t272;
t286 = t213 * t219;
t224 = cos(t228);
t214 = -pkin(1) * cos(t245) - pkin(2) * cos(t233) - pkin(3) * t224;
t285 = t214 * t217;
t225 = cos(t229);
t215 = -pkin(1) * cos(t246) - pkin(2) * cos(t234) - pkin(3) * t225;
t284 = t215 * t218;
t226 = cos(t230);
t216 = -pkin(1) * cos(t247) - pkin(2) * cos(t235) - pkin(3) * t226;
t283 = t216 * t219;
t282 = t217 * t221;
t281 = t217 * t224;
t280 = t218 * t222;
t279 = t218 * t225;
t278 = t219 * t223;
t277 = t219 * t226;
t276 = t231 * t263;
t264 = pkin(1) ^ 2;
t265 = Icges(1,3) + Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t264 + t275) * m(3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2 + t264) * m(2) + 0.2e1 * ((m(2) * rSges(2,1) + m(3) * pkin(2)) * cos(pkin(7)) - m(2) * sin(pkin(7)) * rSges(2,2)) * pkin(1);
t207 = t265 + 0.2e1 * t294;
t206 = t265 + 0.2e1 * t293;
t205 = t265 + 0.2e1 * t292;
t204 = (t210 * t226 + t216 * t276) * t219;
t203 = (t209 * t225 + t215 * t276) * t218;
t202 = (t208 * t224 + t214 * t276) * t217;
t201 = (t210 * t223 + t213 * t276) * t219;
t200 = (t209 * t222 + t212 * t276) * t218;
t199 = (t208 * t221 + t211 * t276) * t217;
t198 = (t207 * t226 + t216 * t289) * t219;
t197 = (t206 * t225 + t215 * t290) * t218;
t196 = (t205 * t224 + t214 * t291) * t217;
t195 = (t207 * t223 + t213 * t289) * t219;
t194 = (t206 * t222 + t212 * t290) * t218;
t193 = (t205 * t221 + t211 * t291) * t217;
t1 = [t196 * t281 + t197 * t279 + t198 * t277 + m(4) + (t202 * t285 + t203 * t284 + t204 * t283) * t263, t196 * t282 + t197 * t280 + t198 * t278 + (t202 * t288 + t203 * t287 + t204 * t286) * t263, 0; t193 * t281 + t194 * t279 + t195 * t277 + (t199 * t285 + t200 * t284 + t201 * t283) * t263, t193 * t282 + t194 * t280 + t195 * t278 + m(4) + (t199 * t288 + t200 * t287 + t201 * t286) * t263, 0; 0, 0, 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
