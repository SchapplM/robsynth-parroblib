% Calculate Gravitation load for parallel robot
% P3PRRR2G2P3A0
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
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR2G2P3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:32
% EndTime: 2020-03-09 21:21:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (276->73), mult. (468->121), div. (45->5), fcn. (285->24), ass. (0->66)
t577 = m(3) / pkin(2);
t542 = sin(qJ(2,3));
t576 = t542 * pkin(1);
t544 = sin(qJ(2,2));
t575 = t544 * pkin(1);
t546 = sin(qJ(2,1));
t574 = t546 * pkin(1);
t538 = legFrame(3,2);
t524 = sin(t538);
t527 = cos(t538);
t512 = t524 * g(1) + t527 * g(2);
t535 = qJ(2,3) + qJ(3,3);
t518 = sin(t535);
t553 = rSges(3,2) * g(3);
t555 = rSges(3,1) * g(3);
t497 = (rSges(3,1) * t512 - t553) * t518 + (rSges(3,2) * t512 + t555) * cos(t535);
t541 = sin(qJ(3,3));
t532 = 0.1e1 / t541;
t573 = t497 * t532;
t539 = legFrame(2,2);
t525 = sin(t539);
t528 = cos(t539);
t513 = t525 * g(1) + t528 * g(2);
t536 = qJ(2,2) + qJ(3,2);
t519 = sin(t536);
t498 = (rSges(3,1) * t513 - t553) * t519 + (rSges(3,2) * t513 + t555) * cos(t536);
t543 = sin(qJ(3,2));
t533 = 0.1e1 / t543;
t572 = t498 * t533;
t540 = legFrame(1,2);
t526 = sin(t540);
t529 = cos(t540);
t514 = t526 * g(1) + t529 * g(2);
t537 = qJ(2,1) + qJ(3,1);
t520 = sin(t537);
t499 = (rSges(3,1) * t514 - t553) * t520 + (rSges(3,2) * t514 + t555) * cos(t537);
t545 = sin(qJ(3,1));
t534 = 0.1e1 / t545;
t571 = t499 * t534;
t531 = g(3) * pkin(1) * m(3);
t548 = cos(qJ(2,3));
t554 = rSges(2,2) * g(3);
t556 = rSges(2,1) * g(3);
t570 = t532 * ((t531 + m(2) * (rSges(2,2) * t512 + t556)) * t548 - m(2) * (-rSges(2,1) * t512 + t554) * t542 + (t512 * t576 + t497) * m(3));
t550 = cos(qJ(2,2));
t569 = t533 * ((t531 + m(2) * (rSges(2,2) * t513 + t556)) * t550 - m(2) * (-rSges(2,1) * t513 + t554) * t544 + (t513 * t575 + t498) * m(3));
t552 = cos(qJ(2,1));
t568 = t534 * ((t531 + m(2) * (rSges(2,2) * t514 + t556)) * t552 - m(2) * (-rSges(2,1) * t514 + t554) * t546 + (t514 * t574 + t499) * m(3));
t567 = t542 * t541;
t566 = t544 * t543;
t565 = t546 * t545;
t547 = cos(qJ(3,3));
t564 = (t548 * t547 - t567) * t570;
t549 = cos(qJ(3,2));
t563 = (t550 * t549 - t566) * t569;
t551 = cos(qJ(3,1));
t562 = (t552 * t551 - t565) * t568;
t561 = ((pkin(2) * t547 + pkin(1)) * t548 - pkin(2) * t567) * t573;
t560 = ((pkin(2) * t549 + pkin(1)) * t550 - pkin(2) * t566) * t572;
t559 = ((pkin(2) * t551 + pkin(1)) * t552 - pkin(2) * t565) * t571;
t558 = 0.1e1 / pkin(1);
t530 = m(1) + m(2) + m(3);
t517 = t529 * g(1) - t526 * g(2);
t516 = t528 * g(1) - t525 * g(2);
t515 = t527 * g(1) - t524 * g(2);
t1 = [-m(4) * g(1) + (-t515 * t527 - t516 * t528 - t517 * t529) * t530 + (t524 * t564 + t525 * t563 + t526 * t562 + (-t524 * t561 - t525 * t560 - t526 * t559) * t577) * t558; -m(4) * g(2) + (t515 * t524 + t516 * t525 + t517 * t526) * t530 + (t527 * t564 + t528 * t563 + t529 * t562 + (-t527 * t561 - t528 * t560 - t529 * t559) * t577) * t558; -m(4) * g(3) + (-t518 * t570 - t519 * t569 - t520 * t568 + ((pkin(2) * t520 + t574) * t571 + (pkin(2) * t519 + t575) * t572 + (pkin(2) * t518 + t576) * t573) * t577) * t558;];
taugX  = t1;
