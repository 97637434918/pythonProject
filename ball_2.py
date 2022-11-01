import pygame
import random
import numpy as np

SCREEN_WIDTH =500
SCREEN_HEIGHT = 400  
FPS = 50
       
class circle(pygame.sprite.Sprite):        
    def __init__(self, screen,r=20):  
        pygame.sprite.Sprite.__init__(self)  
        self.r=r
        self.screen=screen
        self.color=(random.randrange(50,255),random.randrange(50,255),random.randrange(50,255))
        self.image = pygame.Surface((2*r,2*r))
        self.image.set_colorkey((20,120,20))
        self.image.fill((20,120,20))
        self.rect=self.image.get_rect()
        self.rect.x = random.randrange(SCREEN_WIDTH-self.rect.width) 
        self.rect.y = random.randrange(SCREEN_HEIGHT-self.rect.height)
        pygame.draw.circle(self.image,self.color,(self.r,self.r),self.r) 
        self.radius=r  #this makes ".collide_circle" to detect using circle,rather than rect
        self.speedx = random.randrange(1,5)
        self.speedy = random.randrange(1,5) 
    def update(self):
        self.rect.move_ip(self.speedx, self.speedy)
        self.v=np.array([self.speedx,self.speedy])
        if (self.rect.y<=0 or self.rect.y>=SCREEN_HEIGHT-self.rect.height):
            self.speedy=-self.speedy
        if self.rect.x<=0 or self.rect.x>=SCREEN_WIDTH-self.rect.width:
            self.speedx=-self.speedx
        self.screen.blit(self.image, self.rect)
            
     
pygame.init()
clock = pygame.time.Clock()
screen = pygame.display.set_mode((SCREEN_WIDTH , SCREEN_HEIGHT))
c1=circle(screen,30)
balls=pygame.sprite.Group()
for i in range(10):
    c=circle(screen)
    balls.add(c)
run=True
while pygame.sprite.spritecollideany(c1,balls): #避免一开始就碰撞在一起
    c1=circle(screen,30)
while run:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run=False
    c=pygame.sprite.spritecollideany(c1,balls)
    if c:    
        #exchange the speed of c1, and c2
        c1.speedy,c.speedy=c.speedy,c1.speedy
        c1.speedx,c.speedx=c.speedx,c1.speedx
        screen.fill((220, 0, 220))                   
    else:
        screen.fill((0,0,0))
    c1.update()
    balls.update()
    pygame.display.update()
    clock.tick(FPS)
