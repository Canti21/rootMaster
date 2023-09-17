const burguerMenu = document.querySelector('.burguer-menu');
const mobileMenu = document.querySelector('.mobile-menu');
const darken = document.querySelector('.darken')

burguerMenu.addEventListener('click', toggleMobileMenu);

function toggleMobileMenu()
{
    mobileMenu.classList.toggle('inactive')
    darken.classList.toggle('inactive')
}
