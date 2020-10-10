% Calculate Gravitation load for parallel robot
% P3PRRR2G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR2G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:02
% EndTime: 2020-03-09 21:20:02
% DurationCPUTime: 0.21s
% Computational Cost: add. (231->63), mult. (359->103), div. (30->5), fcn. (234->33), ass. (0->46)
t425 = legFrame(1,2);
t413 = sin(t425);
t416 = cos(t425);
t389 = t413 * g(1) + t416 * g(2);
t392 = t416 * g(1) - t413 * g(2);
t422 = qJ(2,1) + qJ(3,1);
t444 = (rSges(3,1) * t392 + rSges(3,2) * t389) * cos(t422) + (rSges(3,1) * t389 - rSges(3,2) * t392) * sin(t422);
t424 = legFrame(2,2);
t412 = sin(t424);
t415 = cos(t424);
t388 = t412 * g(1) + t415 * g(2);
t391 = t415 * g(1) - t412 * g(2);
t421 = qJ(2,2) + qJ(3,2);
t443 = (rSges(3,1) * t391 + rSges(3,2) * t388) * cos(t421) + (rSges(3,1) * t388 - rSges(3,2) * t391) * sin(t421);
t423 = legFrame(3,2);
t411 = sin(t423);
t414 = cos(t423);
t387 = t411 * g(1) + t414 * g(2);
t390 = t414 * g(1) - t411 * g(2);
t420 = qJ(2,3) + qJ(3,3);
t442 = (rSges(3,1) * t390 + rSges(3,2) * t387) * cos(t420) + (rSges(3,1) * t387 - rSges(3,2) * t390) * sin(t420);
t441 = pkin(1) * m(3);
t440 = m(3) / pkin(2);
t417 = 0.1e1 / sin(qJ(3,3));
t439 = t442 * t417;
t418 = 0.1e1 / sin(qJ(3,2));
t438 = t443 * t418;
t419 = 0.1e1 / sin(qJ(3,1));
t437 = t444 * t419;
t430 = t417 * ((t390 * t441 + m(2) * (rSges(2,1) * t390 + rSges(2,2) * t387)) * cos(qJ(2,3)) - (-t387 * t441 + m(2) * (-rSges(2,1) * t387 + rSges(2,2) * t390)) * sin(qJ(2,3)) + t442 * m(3));
t429 = t418 * ((t391 * t441 + m(2) * (rSges(2,1) * t391 + rSges(2,2) * t388)) * cos(qJ(2,2)) - (-t388 * t441 + m(2) * (-rSges(2,1) * t388 + rSges(2,2) * t391)) * sin(qJ(2,2)) + t443 * m(3));
t428 = t419 * ((t392 * t441 + m(2) * (rSges(2,1) * t392 + rSges(2,2) * t389)) * cos(qJ(2,1)) - (-t389 * t441 + m(2) * (-rSges(2,1) * t389 + rSges(2,2) * t392)) * sin(qJ(2,1)) + t444 * m(3));
t402 = -t423 + qJ(2,3);
t403 = -t424 + qJ(2,2);
t404 = -t425 + qJ(2,1);
t427 = 0.1e1 / pkin(1);
t401 = qJ(3,1) + t404;
t400 = qJ(3,2) + t403;
t399 = qJ(3,3) + t402;
t398 = cos(t401);
t397 = cos(t400);
t396 = cos(t399);
t395 = sin(t401);
t394 = sin(t400);
t393 = sin(t399);
t1 = [-m(4) * g(1) + (-t393 * t430 - t394 * t429 - t395 * t428 + ((pkin(2) * t395 + pkin(1) * sin(t404)) * t437 + (pkin(2) * t394 + pkin(1) * sin(t403)) * t438 + (pkin(2) * t393 + pkin(1) * sin(t402)) * t439) * t440) * t427; -m(4) * g(2) + (t396 * t430 + t397 * t429 + t398 * t428 + ((-pkin(2) * t398 - pkin(1) * cos(t404)) * t437 + (-pkin(2) * t397 - pkin(1) * cos(t403)) * t438 + (-pkin(2) * t396 - pkin(1) * cos(t402)) * t439) * t440) * t427; -0.3e1 * g(3) * (m(1) + m(2) + m(3)) - m(4) * g(3);];
taugX  = t1;
