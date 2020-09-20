% Calculate Gravitation load for parallel robot
% P3PRRRR1G3A0
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
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:18
% EndTime: 2020-03-09 21:02:18
% DurationCPUTime: 0.39s
% Computational Cost: add. (252->70), mult. (525->141), div. (45->10), fcn. (366->18), ass. (0->56)
t552 = sin(qJ(3,1));
t558 = cos(qJ(3,1));
t583 = -rSges(3,1) * t558 + rSges(3,2) * t552;
t550 = sin(qJ(3,2));
t556 = cos(qJ(3,2));
t582 = -rSges(3,1) * t556 + rSges(3,2) * t550;
t548 = sin(qJ(3,3));
t554 = cos(qJ(3,3));
t581 = -rSges(3,1) * t554 + rSges(3,2) * t548;
t580 = rSges(3,1) * g(3);
t545 = legFrame(3,2);
t532 = sin(t545);
t535 = cos(t545);
t526 = t532 * g(1) + t535 * g(2);
t549 = sin(qJ(2,3));
t539 = 0.1e1 / t549;
t573 = t526 * t539;
t546 = legFrame(2,2);
t533 = sin(t546);
t536 = cos(t546);
t527 = t533 * g(1) + t536 * g(2);
t551 = sin(qJ(2,2));
t540 = 0.1e1 / t551;
t572 = t527 * t540;
t547 = legFrame(1,2);
t534 = sin(t547);
t537 = cos(t547);
t528 = t534 * g(1) + t537 * g(2);
t553 = sin(qJ(2,1));
t541 = 0.1e1 / t553;
t571 = t528 * t541;
t570 = t539 * t548;
t569 = t540 * t550;
t568 = t541 * t552;
t529 = t535 * g(1) - t532 * g(2);
t555 = cos(qJ(2,3));
t523 = ((-rSges(2,1) * t526 + rSges(2,2) * t529) * t555 + t549 * (rSges(2,1) * t529 + rSges(2,2) * t526)) * m(2) + ((-t529 * rSges(3,3) + t581 * t526) * t555 + t549 * (-t526 * rSges(3,3) - t581 * t529)) * m(3);
t542 = 0.1e1 / t554;
t567 = t523 * t539 * t542;
t530 = t536 * g(1) - t533 * g(2);
t557 = cos(qJ(2,2));
t524 = ((-rSges(2,1) * t527 + rSges(2,2) * t530) * t557 + t551 * (rSges(2,1) * t530 + rSges(2,2) * t527)) * m(2) + ((-t530 * rSges(3,3) + t582 * t527) * t557 + t551 * (-t527 * rSges(3,3) - t582 * t530)) * m(3);
t543 = 0.1e1 / t556;
t566 = t524 * t540 * t543;
t531 = t537 * g(1) - t534 * g(2);
t559 = cos(qJ(2,1));
t525 = ((-rSges(2,1) * t528 + rSges(2,2) * t531) * t559 + t553 * (rSges(2,1) * t531 + rSges(2,2) * t528)) * m(2) + ((-t531 * rSges(3,3) + t583 * t528) * t559 + t553 * (-t528 * rSges(3,3) - t583 * t531)) * m(3);
t544 = 0.1e1 / t558;
t565 = t525 * t541 * t544;
t564 = t549 * t526 + t529 * t555;
t563 = t551 * t527 + t530 * t557;
t562 = t553 * t528 + t531 * t559;
t561 = 0.1e1 / pkin(2);
t560 = rSges(3,2) * g(3);
t538 = m(1) + m(2) + m(3);
t1 = [-m(4) * g(1) + (-t535 * t567 - t536 * t566 - t537 * t565) * t561 + (-(t553 * t534 + t537 * t559) * t571 - (t551 * t533 + t536 * t557) * t572 - (t549 * t532 + t535 * t555) * t573) * t538; -m(4) * g(2) + (t532 * t567 + t533 * t566 + t534 * t565) * t561 + (-(-t534 * t559 + t537 * t553) * t571 - (-t533 * t557 + t536 * t551) * t572 - (-t532 * t555 + t535 * t549) * t573) * t538; -m(4) * g(3) + (-t542 * t526 * t570 - t543 * t527 * t569 - t544 * t528 * t568) * t538 + (-t559 / t558 ^ 2 * t525 * t568 - t557 / t556 ^ 2 * t524 * t569 - t555 / t554 ^ 2 * t523 * t570 + (t544 * ((t562 * rSges(3,2) - t580) * t558 + t552 * (t562 * rSges(3,1) + t560)) + t543 * ((t563 * rSges(3,2) - t580) * t556 + t550 * (t563 * rSges(3,1) + t560)) + t542 * ((t564 * rSges(3,2) - t580) * t554 + t548 * (t564 * rSges(3,1) + t560))) * m(3)) * t561;];
taugX  = t1;
